function [NewPopulation,NewFitness] = AngleConstrainedSelection(Population,N,Z,W)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yubo Wang (email: 1202111230@cug.edu.cn)

    PopObj  = Population.objs - Z;
    Nonzero = any(PopObj ~= 0,2);
    % A solution at the ideal point has equal priority in every direction.
    Angle  = zeros(size(PopObj,1),size(W,1));
    Cosine = 1-pdist2(PopObj(Nonzero,:),W,'cosine');
    Angle(Nonzero,:) = acos(max(-1,min(1,Cosine)));
    [~,Associate]    = min(Angle,[],2);
    NewPopulation    = [];
    NewFitness       = [];
    Fitness          = CalFitness(Population.objs,Population.cons);
    for i = 1 : N
        Candidate    = Population(Associate == i);
        CandidateFit = Fitness(Associate == i);
        if isempty(Candidate)
            Theta = rad2deg(Angle(:,i));
            Index = find(Theta==min(Theta));
            Index = Index(1);
            NewPopulation = [NewPopulation Population(Index)];
            NewFitness   = [NewFitness Fitness(Index)];
        else
            CV = sum(max(0,Candidate.cons),2);
            if length(Candidate) == 1
                NewPopulation = [NewPopulation Candidate];
                NewFitness    = [NewFitness CandidateFit];
            else
                FeasiblePop    = Candidate(CV ==0);
                FeasiblePopFit = CandidateFit(CV ==0);
                if length(FeasiblePop) == 1
                    NewPopulation = [NewPopulation FeasiblePop];
                    NewFitness    = [NewFitness FeasiblePopFit];
                elseif length(FeasiblePop) > 1
                    Index = find(FeasiblePopFit==min(FeasiblePopFit));
                    Index = Index(1);
                    NewPopulation = [NewPopulation FeasiblePop(Index)];
                    NewFitness    = [NewFitness FeasiblePopFit(Index)];
                else
                    Index = find(CV==min(CV));
                    Index = Index(1);
                    NewPopulation = [NewPopulation Candidate(Index)];
                    NewFitness    = [NewFitness CandidateFit(Index)];
                end
            end
        end
    end
end