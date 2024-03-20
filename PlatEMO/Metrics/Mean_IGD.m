function score = Mean_IGD(Population,Problem)
% <min> <multi/many> <real/integer/label/binary/permutation> <large/none> <constrained/none> <expensive/none> <multimodal/none> <sparse/none> <dynamic/none> <robust>
% Mean inverted generational distance (for robust optimization)

%------------------------------- Reference --------------------------------
% C. A. Coello Coello and N. C. Cortes, Solving multiobjective optimization
% problems using an artificial immune system, Genetic Programming and
% Evolvable Machines, 2005, 6(2): 163-190.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N     = length(Population);
    PopX  = Problem.Perturb(Population.decs);
    score = [];
    for i = 1 : N : length(PopX)
        PopObj = PopX(i:i+N-1).best.objs;
        if size(PopObj,2) ~= size(Problem.optimum,2)
            score = [score,nan];
        else
            score = [score,mean(min(pdist2(Problem.optimum,PopObj),[],2))];
        end
    end
    score = mean(score(~isnan(score)));
end