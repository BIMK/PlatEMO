function [ Offspring ] = OperatorDE_pbest_1_main(Population, popsize, Problem, Fitness, p)
% The operator of DE

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Kangjia Qiao (email: qiaokangjia@yeah.net)

    permutation = randperm(Problem.N);
    r0 = permutation;
    [r1,r2,r3] = gnR1R2R3(Problem.N, r0);

    array = permutation(1:popsize);
    pop1 = Population(array);

    [~, indBest] = sort(Fitness, 'ascend');
    pNP = max(round(p * Problem.N), 2);         % choose at least two best solutions  
    randindex = ceil(rand(1, popsize) * pNP);	% select from [1, 2, 3, ..., pNP]
    randindex = max(1, randindex);              % to avoid the problem that rand = 0 and thus ceil(rand) = 0
    pbest = Population(indBest(randindex));     % randomly choose one of the top 100p% solutions

    Offspring = OperatorDE_pbest_1(Problem,Population(array),pbest,Population(r1(1:popsize)),Population(r2(1:popsize)));
end