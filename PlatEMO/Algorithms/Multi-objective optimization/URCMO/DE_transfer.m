function [ offspring ] = DE_transfer(Problem, Population1, Population2, popsize)

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

    index = randi([1,length(Fm)],popsize,1);
    F     = Fm(index);
    F     = F';
    index = randi([1,length(CRm)],popsize,1);
    CR    = CRm(index);
    CR    = CR';

    index =randi(Problem.N, popsize,1);

    permutation = randperm(Problem.N);

    array = permutation(1:popsize);

    pop1 = Population1(array).decs;
    pop2 = Population2.decs;

    vi = pop2(index,:);


    mask  = rand(popsize, Problem.D) > CR(:, ones(1, Problem.D)); % mask is used to indicate which elements of ui comes from the parent
    rows  = (1 : popsize)'; cols = floor(rand(popsize, 1) * Problem.D)+1; % choose one position where the element of ui doesn't come from the parent
    jrand = sub2ind([popsize Problem.D], rows, cols); mask(jrand) = false;
    u     = vi; u(mask) = pop1(mask);

    offspring = Problem.Evaluation(u);
end