function [ Offspring ] = OperatorDE_pbest_1_main_OP(Population, popsize, Problem, Fitness, DV_OP, DV_NOP, p)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Kangjia Qiao (email: qiaokangjia@yeah.net)

    permutation = randperm(Problem.N);
    r0          = permutation;
    [r1,r2,~]   = gnR1R2R3(Problem.N, r0);
    
    array  = permutation(1:popsize);
    OffDec = Population(array).decs;
    
    [~,indBest] = sort(Fitness, 'ascend');
    pNP         = max(round(p * Problem.N), 2);     % choose at least two best solutions
    randindex   = ceil(rand(1, popsize) * pNP);     % select from [1, 2, 3, ..., pNP]
    randindex   = max(1, randindex);                % to avoid the problem that rand = 0 and thus ceil(rand) = 0
    pbest       = Population(indBest(randindex));   % randomly choose one of the top 100p% solutions
    
    NewDec = OperatorDE_pbest_1(Problem,Population(array).decs,pbest.decs,Population(r1(1:popsize)).decs,Population(r2(1:popsize)).decs);
    OffDec(:,DV_OP) = NewDec(:,DV_OP);
    
    %% DV_NOP: Randomly select individuals to assign dimensional information
    if sum(DV_NOP)>1
        rand_ind = randi([1 popsize],popsize,length(DV_NOP));
        Decc     = Population.decs;
        for h = 1 : popsize
            Ori = Decc(rand_ind(h,:), DV_NOP);
            Oth(h,:) = diag(Ori);
        end
        OffDec(:,DV_NOP) = Oth;
    end
    Offspring = OffDec;
    Offspring = Problem.Evaluation(Offspring);
end