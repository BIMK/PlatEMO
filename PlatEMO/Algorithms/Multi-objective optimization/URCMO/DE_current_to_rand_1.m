function [ offspring ] = DE_current_to_rand_1(Problem, Population, popsize)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Kangjia Qiao

    Fm  = [0.6,0.8,1.0];
    CRm = [0.1,0.2,1.0];
    lu  = [Problem.lower;Problem.upper];

    index = randi([1,length(Fm)],popsize,1);
    F     = Fm(index);
    F     = F';
    index = randi([1,length(CRm)],popsize,1);
    
    permutation = randperm(Problem.N);
    r0          = permutation;
    [r1, r2,r3] = gnR1R2R3(Problem.N,  r0);

    array = permutation(1:popsize);
    pop   = Population.decs;

    pop1 = pop(1:popsize,:);

    vi = pop1 +  repmat(rand(popsize,1),1,Problem.D) .*(pop(r1(array),:) - pop1 ) + F(:, ones(1, Problem.D)) .*(pop(r2(array),:) - pop(r3(array),:));
    vi = boundConstraint(vi, pop1, lu);
    u  = vi;
    
    offspring = Problem.Evaluation(u);
end