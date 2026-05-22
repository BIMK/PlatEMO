function Reward = EstimateReward(Problem, ArchiveOld, ArchiveNew, ratio,Cons0,score_HV0,optimum0)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Phi     = Problem.FE/Problem.maxFE;
    PopObj1 = best(ArchiveOld);
    PopObj2 = best(ArchiveNew);

    if isempty(PopObj1) || isempty(PopObj2)
        Cons1  = sum(sum(max(ArchiveOld.cons,0),2));
        Cons2  = sum(sum(max(ArchiveNew.cons,0),2));
        Reward = (Cons1-Cons2)/Cons1;
    else
        PopObj1 = ArchiveOld(all(ArchiveOld.cons<=0,2));
        PopObj2 = ArchiveNew(all(ArchiveNew.cons<=0,2));

        [~,k] = size(PopObj1);
        for i = 1 : k
            obj1(i,:) = PopObj1(i).obj;
        end
        [~,k] = size(PopObj2);
        for i = 1 : k
            obj2(i,:) = PopObj2(i).obj;
        end
        [~,m] = size(obj1);
        for i = 1 : m
            optimum1(i) = max(obj1(:,i));
            optimum2(i) = max(obj2(:,i));
        end
        optimum1  = optimum1 * 1.1;
        score_HV1 = HV(ArchiveOld,optimum1);
        score_HV2 = HV(ArchiveNew,optimum1);
        score_HV  = (score_HV2 - score_HV1)/(score_HV1);

        score_PD1 = PD(ArchiveOld);
        score_PD2 = PD(ArchiveNew);
        score_PD  = (score_PD2 - score_PD1)/(score_PD1);

        Reward = Phi * score_PD + (1-Phi) * score_HV;
    end
end