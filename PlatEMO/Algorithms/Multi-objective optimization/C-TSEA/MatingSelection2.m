function MatingPool = MatingSelection2(Population,Archive,N)
% The mating selection of stage 2

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming

    Fitness1    = CalFitness(Population.objs,Population.cons);
    MatingPool1 = TournamentSelection(2,N,Fitness1);
    MatingPool  = Population(MatingPool1);
end