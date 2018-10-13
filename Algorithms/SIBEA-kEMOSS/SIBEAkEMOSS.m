function SIBEAkEMOSS(Global)
% <algorithm> <O-Z>
% Improving Hypervolume-based Multiobjective Evolutionary Algorithms by
% Using Objective Reduction Methods
% G --- 5 --- Reduction frequency generations
% k --- 2 --- Size of reduced objective set

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
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
        Offspring     = Global.Variation(Population(MatingPool));
        Population    = EnvironmentalSelection([Population,Offspring],Global.N,objective_set);
    end
end