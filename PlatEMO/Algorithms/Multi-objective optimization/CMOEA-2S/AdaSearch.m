function Offspring = AdaSearch(Problem, Population, Fitness, IsTwo)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    if IsTwo
        MatingPool  = TournamentSelection(2,Problem.N,Fitness);
        MatingPool0 = TournamentSelection(2,Problem.N/2,Fitness);
        Offspring   = OperatorDE(Problem,Population(MatingPool0),Population(MatingPool(1:end/2)),Population(MatingPool(end/2+1:end)));
    else
        MatingPool1 = TournamentSelection(2,Problem.N/2,Fitness);
        MatingPool2 = TournamentSelection(2,Problem.N/2,Fitness);
        MatingPool0 = TournamentSelection(2,Problem.N/2,Fitness);
        Offspring   = OperatorDE(Problem,Population(MatingPool0),Population(MatingPool1),Population(MatingPool2));
    end
end