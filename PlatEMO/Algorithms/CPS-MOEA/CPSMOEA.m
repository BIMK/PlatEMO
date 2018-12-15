function CPSMOEA(Global)
% <algorithm> <C>
% Classification and Pareto domination based multi-objective evolutionary
% algorithm
% M --- 3 --- Number of generated offsprings for each solution

%------------------------------- Reference --------------------------------
% J. Zhang, A. Zhou, and G. Zhang, A classification and Pareto domination
% based multiobjective evolutionary algorithm, Proceedings of the IEEE
% Congress on Evolutionary Computation, 2015, 2883-2890.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    M = Global.ParameterSet(3);

    %% Generate random population
    Population   = Global.Initialization();
    [Pgood,Pbad] = NDS(Population,floor(Global.N/2));

    %% Optimization
    while Global.NotTermination(Population)
        KNN(Pgood.decs,Pbad.decs);
        Offspring  = Operator(Population,M);
        Population = NDS([Population,Offspring],Global.N);
        FrontNo    = NDSort(Offspring.objs,1);
        Pgood      = NDS([Pgood,Offspring(FrontNo==1)],floor(Global.N/2));
        Pbad       = NDS([Pbad,Offspring(FrontNo~=1)],floor(Global.N/2));
    end
end