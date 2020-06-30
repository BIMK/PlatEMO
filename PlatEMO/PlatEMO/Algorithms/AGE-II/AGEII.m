function AGEII(Global)
% <algorithm> <A>
% Approximation-guided evolutionary multi-objective algorithm II
% epsilon --- 0.1 --- The parameter in grid location calculation

%------------------------------- Reference --------------------------------
% M. Wagner and F. Neumann, A fast approximation-guided evolutionary
% multi-objective algorithm, Proceedings of the 15th Annual Conference on
% Genetic and Evolutionary Computation, 2013, 687-694.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    epsilon = Global.ParameterSet(0.1);

    %% Generate the sampling points and random population
    Population = Global.Initialization();
    Archive    = UpdateArchive(Population,epsilon);

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = MatingSelection(Population.objs);
        Offspring  = GA(Population(MatingPool));
        Archive    = UpdateArchive([Archive,Offspring],epsilon);
        Population = EnvironmentalSelection([Population,Offspring],Archive.objs,Global.N);
    end
end