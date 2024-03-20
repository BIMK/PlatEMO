function score = Task1_HV(Population,optimum)
% <max> <multi/many> <real/integer/label/binary/permutation> <large/none> <constrained/none> <expensive/none> <sparse/none> <multitask>
% Hypervolume of the first task (for multitask optimization)

%------------------------------- Reference --------------------------------
% E. Zitzler and L. Thiele, Multiobjective evolutionary algorithms: A
% comparative case study and the strength Pareto approach, IEEE
% Transactions on Evolutionary Computation, 1999, 3(4): 257-271.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    PopDec = Population.decs;
    score  = HV(Population(PopDec(:,end)==1),optimum);
end