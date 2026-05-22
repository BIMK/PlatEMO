function [return_pop,return_Fitness,SelInd] = DI_Update(MaxPop,N,VAR1,VAR2,W)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    PopObj  = MaxPop.objs;
    Zmin    = min(PopObj,[],1);
    [Num,~] = size(PopObj); 
    Cons    = max(0,MaxPop.cons);
    CV      = sum(Cons,2);
    Gmin    = min(CV);
    Gmax    = max(CV);
    
    %% shift the objective space to R+
    PopObj = PopObj - repmat(Zmin,Num,1) + 1e-6;
    
    findex       = find(CV<=VAR1);
    ifindex      = find(CV>VAR1);
    fPopulation  = MaxPop(findex);
    ifPopulation = MaxPop(ifindex);
    Population   = [fPopulation,ifPopulation];  
    fnum         = length(fPopulation);
    inum         = length(ifPopulation);
    fPopObj      = PopObj(findex,:);
    ifPopObj     = PopObj(ifindex,:);
    PopObj       = [fPopObj;ifPopObj]; 
    fCV          = CV(findex);
    iCV          = CV(ifindex);
    CV           = [fCV;iCV];
    
    %% calculate the indicator matrix
    IMatrix = ones(Num,Num);
    for i = 1 : Num
        Ci = CV(i);
        if Ci <= VAR2
            Fi           = PopObj(i,:);
            Ir           = log(repmat(Fi,fnum,1)./fPopObj);
            MaxIr        = max(Ir,[],2);
            MinIr        = min(Ir,[],2);
            CVA          = MaxIr;
            DomInds      = find(MaxIr<=0);
            CVA(DomInds) = MinIr(DomInds);
            IndicatorV   = CVA;
        else
            IC         = repmat(Ci+1e-6,fnum,1)./(fCV+1e-6); 
            Fi         = PopObj(i,:);
            MaxF       = max(repmat(Fi,fnum,1),fPopObj);
            MinF       = min(repmat(Fi,fnum,1),fPopObj);
            CVF        = max(MaxF./MinF,[],2);
            IndicatorV = log(max([CVF,IC],[],2));
        end
        IMatrix(1:fnum,i) = IndicatorV;
        
        
        IC                    = repmat(Ci+1e-6,inum,1)./(iCV+1e-6); 
        IndicatorV            = log(IC)+repmat(log((Gmin+1e-6)/(Gmax+1e-6))+log(min(min(PopObj,[],1)./max(PopObj,[],1))),inum,1)-1e-6;
        IMatrix(1+fnum:Num,i) = IndicatorV;
    end
    
    IMatrix(logical(eye(Num))) = Inf;
    Fitness                    = min(IMatrix,[],2);

    % using indicator-based CHT to update the population
    SelInd         = Indicator_based_CHT(PopObj,IMatrix,W,N);
    return_pop     = Population(SelInd);
    return_Fitness = Fitness(SelInd);
end