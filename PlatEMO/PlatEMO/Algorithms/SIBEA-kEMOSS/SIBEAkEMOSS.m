function SIBEAkEMOSS(Global)
% <algorithm> <S>
% SIBEA with minimum objective subset of size k with minimum error
% G --- 5 --- Reduction frequency generations
% k --- 2 --- Size of reduced objective set

%------------------------------- Reference --------------------------------
% D. Brockhoff and E. Zitzler, Improving hypervolume-based multiobjective
% evolutionary algorithms by using objective reduction methods, Proceedings
% of the IEEE Congress on Evolutionary Computation, 2007, 2086-2093.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Liangli Zhen

   %% Parameter setting
   [G,k] = Global.ParameterSet(5,2);
    
    %% Generate random population
    Population    = Global.Initialization();
    iteration_num = 0;
    objective_set = 1 : k;

    %% Optimization
    while Global.NotTermination(Population)
        if mod(iteration_num,G) ==0
            objective_set = kEMOSS(Population,k);
        end
        iteration_num = iteration_num + 1;
        MatingPool    = randi(length(Population),1,Global.N);
        Offspring     = GA(Population(MatingPool));
        Population    = EnvironmentalSelection([Population,Offspring],Global.N,objective_set);
    end
end