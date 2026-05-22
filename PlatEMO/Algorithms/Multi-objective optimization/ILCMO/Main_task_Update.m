function [return_pop,return_Fitness,SelInd] = Main_task_Update(MaxPop,N,Fmin)

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
    % shift the objective space to R+
    PopObj = PopObj - repmat(Zmin,Num,1) + 1e-6;
    
    findex       = find(CV<=0);
    ifindex      = find(CV>0);
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
        if Ci <= 0 %%%%% Xi is feasible
            Fi           = PopObj(i,:);
            Ir           = log(repmat(Fi,fnum,1)./fPopObj);
            MaxIr        = max(Ir,[],2);
            MinIr        = min(Ir,[],2);
            CVA          = MaxIr;
            DomInds      = find(MaxIr<=0);
            CVA(DomInds) = MinIr(DomInds);
            IndicatorV   = CVA;
        else  %%%%% Xi is an infeasible solution
            IC         = Inf(fnum,1); 
            IndicatorV = IC;
        end
        IMatrix(1:fnum,i)     = IndicatorV;
        IC                    = repmat(Ci+1e-6,inum,1)./(iCV+1e-6); 
        IndicatorV            = log(IC)+repmat(log((Gmin+1e-6)/(Gmax+1e-6))+log(min(min(PopObj,[],1)./max(PopObj,[],1))),inum,1)-1e-6;
        IMatrix(1+fnum:Num,i) = IndicatorV;
    end
    
    IMatrix(logical(eye(Num))) = Inf;
    IrFitness                  = min(IMatrix,[],2);
    Level1Index                = find(IrFitness>=0);
    Len_Level1                 = length(Level1Index);

    %% update the Archive
    if Len_Level1 <= N
        [~,SortInd]    = sort(-IrFitness);
        SelInd         = SortInd(1:N);
        return_pop     = Population(SelInd);
        return_Fitness = IrFitness(SelInd);
    else
        PopObj         = PopObj(Level1Index,:) + repmat(Zmin,Len_Level1,1) - repmat(Fmin,Len_Level1,1);
        SelInd         = Selection_Operator_of_PREA(PopObj,IMatrix(Level1Index,Level1Index),N);
        return_pop     = Population(Level1Index(SelInd));
        return_Fitness = IrFitness(Level1Index(SelInd));
    end
end