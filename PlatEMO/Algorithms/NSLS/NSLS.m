function NSLS(Global)
% <algorithm> <N>
% Multiobjective optimization framework based on nondominated sorting and
% local search

%------------------------------- Reference --------------------------------
% B. Chen, W. Zeng, Y. Lin, and D. Zhang, A new local search-based
% multiobjective optimization algorithm, IEEE Transactions on Evolutionary
% Computation, 2015, 19(1): 50-73.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate random population
    Population = Global.Initialization();

    %% Optimization
    while Global.NotTermination(Population)
        Offspring  = Operator(Population,{5555,5555});
        Population = EnvironmentalSelection([Population,Offspring],Global.N);
    end
end