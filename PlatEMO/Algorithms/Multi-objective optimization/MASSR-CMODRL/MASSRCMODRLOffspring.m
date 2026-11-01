function Offspring = MASSRCMODRLOffspring(Problem,Population,Fitness,operator,count)
%MASSRCMODRLOffspring - Evaluate the requested number of offspring.
%
%   operator 1 uses DE; operator 2 uses GA. Numeric decision matrices avoid
%   implicit evaluations, and count limits the explicit evaluations at stage
%   and budget boundaries.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    if count == 0
        Offspring = Population([]);
        return;
    end
    if operator == 1
        parents = TournamentSelection(2,2*count,Fitness);
        target  = Population(randperm(numel(Population),count));
        Dec     = OperatorDE(Problem,target.decs,Population(parents(1:count)).decs, ...
            Population(parents(count+1:end)).decs);
    else
        parents = TournamentSelection(2,2*ceil(count/2),Fitness);
        Dec     = OperatorGA(Problem,Population(parents).decs);
    end
    Offspring = Problem.Evaluation(Dec(1:count,:));
end