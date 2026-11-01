function [Population,Fitness] = NicheSubsetSelection(Population,N,RefV,epsilon)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yubo Wang (email: 1202111230@cug.edu.cn)

    R            = size(RefV, 1);
    numSolutions = length(Population);
    PopObjs      = Population.objs;
    PopCons      = Population.cons;
    Fit          = CalFitness(Population.objs,Population.cons,epsilon);
    Cap          = CalculateNicheCapacity(N, R);
    CV           = sum(max(0, PopCons), 2);
    FeasibleIdx  = find(CV <= epsilon);
    if ~isempty(FeasibleIdx)
        [FrontNo, ~]        = NDSort(PopObjs(FeasibleIdx, :),1);
        NonDominatedIndices = FeasibleIdx(FrontNo == 1);
        SelectedIndices     = NonDominatedIndices; 
        RemainingIndices    = setdiff(1:numSolutions, SelectedIndices);
        [A_associate, ~]    = NicheAssociation(PopObjs(SelectedIndices, :), RefV);
        for r = 1 : R 
            NicheSolutions_A_logical = (A_associate == r);
            NicheSolutions_A_Idx     = SelectedIndices(NicheSolutions_A_logical);
            Count_r = length(NicheSolutions_A_Idx);
            Cap_r   = Cap(r);
            if Count_r > Cap_r 
                NumToMove        = Count_r - Cap_r;
                CandidatePop     = Population(NicheSolutions_A_Idx);
                DensityScores    = CalculateLocalDensity(CandidatePop.objs);
                [~, Rank]        = sort(DensityScores, 'descend');
                LocalIdxToMove   = NicheSolutions_A_Idx(Rank(1:NumToMove));
                SelectedIndices  = setdiff(SelectedIndices, LocalIdxToMove);
                RemainingIndices = [RemainingIndices, LocalIdxToMove'];
                [A_associate, ~] = NicheAssociation(PopObjs(SelectedIndices, :), RefV);
            end
        end
    else
        NonDominatedIndices = [];
        SelectedIndices     = NonDominatedIndices;
        RemainingIndices    = setdiff(1:numSolutions, SelectedIndices);
    end
    SelectedPopulation = Population(SelectedIndices);
    RemainingToFill    = N - length(SelectedPopulation);
    if RemainingToFill > 0
        RemainingPopulation = Population(RemainingIndices);
        RemainingFitness    = Fit(RemainingIndices);
        if ~isempty(SelectedIndices)
            [~, SelectedCounts] = NicheAssociation(PopObjs(SelectedIndices, :), RefV);
        else
            SelectedCounts = zeros(R,1);
        end
        [B_associate, ~] = NicheAssociation(RemainingPopulation.objs, RefV);
        for r = 1 : R 
            if RemainingToFill <= 0
                break;
            end
            Needed = Cap(r) - SelectedCounts(r);
            if Needed <= 0
                continue;
            end
            CandidatesIdx_B = find(B_associate == r);
            Candidates      = RemainingPopulation(CandidatesIdx_B);
            Candidates_Fit  = RemainingFitness(CandidatesIdx_B);
            if Needed > length(Candidates)
                Selected = Candidates;
                Needed   = Needed - length(Candidates);
                RefV_r   = RefV(r, :);
                SelectedPopulation = [SelectedPopulation, Selected];
                AngleDist          = acos(1 - pdist2(RemainingPopulation.objs, RefV_r, 'cosine'));
                [~, RankAngle]     = sort(AngleDist);
                Selected           = RemainingPopulation(RankAngle(1:Needed));
                SelectedPopulation = [SelectedPopulation, Selected];
            else
                [~, RankFit]       = sort(Candidates_Fit); 
                Selected           = Candidates(RankFit(1:Needed));
                SelectedPopulation = [SelectedPopulation, Selected];
            end
        end
    end
    Population = SelectedPopulation(1:N);
    Fitness    = CalFitness(Population.objs, Population.cons, epsilon);
end

function Capacity = CalculateNicheCapacity(N, R)
    capacity_base = floor(N / R);
    capacity_rem  = N - (R - 1) * capacity_base;
    Capacity      = repmat(capacity_base, R, 1);
    Capacity(R)   = capacity_rem;
end

function [associate, NicheCount] = NicheAssociation(PopObj, RefV)
    if isempty(PopObj)
        associate  = [];
        NicheCount = zeros(1, size(RefV, 1));
        return;
    end
    PopObj_Norm    = (PopObj - min(PopObj)) ./ (max(PopObj) - min(PopObj) + 1e-8);
    Cosine         = 1 - pdist2(PopObj_Norm, RefV, 'cosine');
    [~, associate] = max(Cosine, [], 2);
    NicheCount     = histcounts(associate, 1:size(RefV,1)+1);
end

function DensityScores = CalculateLocalDensity(PopObj)
    numSolutions   = size(PopObj, 1);
    Distance_local = pdist2(PopObj, PopObj);
    Distance_local(logical(eye(numSolutions))) = inf;
    Distance_local = sort(Distance_local, 2);
    k = floor(sqrt(numSolutions));
    if k < 1
        k = 1;
    end
    DensityScores = 1./(Distance_local(:, k) + 2);
end