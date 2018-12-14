function PESAII(Global)
% <algorithm> <P>
% Pareto envelope-based selection algorithm II
% div --- 10 --- The number of divisions in each objective

%------------------------------- Reference --------------------------------
% D. W. Corne, N. R. Jerram, J. D. Knowles, and M. J. Oates, PESA-II:
% Region-based selection in evolutionary multiobjective optimization,
% Proceedings of the 3rd Annual Conference on Genetic and Evolutionary
% Computation, 2001, 283-290.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    div = Global.ParameterSet(10);

    %% Generate random population
    Population = Global.Initialization();
    
    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = MatingSelection(Population.objs,Global.N,div);
        Offspring  = GA(Population(MatingPool));
        Population = EnvironmentalSelection([Population,Offspring],Global.N,div);
    end
end