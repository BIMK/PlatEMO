function [Population,Archive] = ICMA_Update(MaxPop,N,W,Zmin,Fmin)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Jiawei Yuan

    PopObj               = MaxPop.objs;
    [Num,M]              = size(PopObj);
    Cons                 = max(0,MaxPop.cons);
    NormCons             = Cons./repmat(max(1,max(Cons,[],1)),Num,1);

    CV                   = sum(NormCons,2);

    % shift the objective space to R+
    PopObj               = PopObj - repmat(Zmin,Num,1) + 1e-6;

    % calculate the indicator matrix
    IMatrix              = ones(Num,Num);
    for i = 1:1:Num
        Ci               = CV(i);    
        if Ci == 0 %%%%% Xi is feasible
            Fi               = PopObj(i,:);
            Ir               = log(repmat(Fi,Num,1)./PopObj);
            MaxIr            = max(Ir,[],2);
            MinIr            = min(Ir,[],2);
            CVA              = MaxIr;
            DomInds          = find(MaxIr<=0);
            CVA(DomInds)     = MinIr(DomInds);
            IndicatorV       = CVA;
        else  %%%%% Xi is an infeasible solution
            IC               = repmat(Ci+1e-6,Num,1)./(CV+1e-6); 
            Fi = PopObj(i,:);
            MaxF = max(repmat(Fi,Num,1),PopObj);
            MinF = min(repmat(Fi,Num,1),PopObj);
            CVF = max(MaxF./MinF,[],2);
            IndicatorV       = log(max([CVF,IC],[],2));
        end
        IMatrix(:,i)     = IndicatorV;
        IMatrix(i,i)     = Inf;
    end

    FeasibleInd                   = find(CV==0);
    Len_F                         = length(FeasibleInd);

    if Len_F<=N
        [~,CV_SortInd]            = sort(CV);
        Archive                   = MaxPop(CV_SortInd(1:N));
    else
        FPopObj                   = PopObj(FeasibleInd,:) + repmat(Zmin,Len_F,1) - repmat(Fmin,Len_F,1);
        SelInd                    = Selection_Operator_of_PREA(FPopObj,IMatrix(FeasibleInd,FeasibleInd),N);
        Archive                   = MaxPop(FeasibleInd(SelInd));
    end

    %%%%%%%%%%%%%  using indicator-based CHT to update the population
    SelInd                        = Indicator_based_CHT(PopObj,IMatrix,W,N);
    Population                    = MaxPop(SelInd);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end