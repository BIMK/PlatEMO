function CCMO(Global)
% <algorithm> <C>
% Coevolutionary constrained multi-objective optimization framework
% type --- 1 --- Type of operator (1. GA 2. DE)

%------------------------------- Reference --------------------------------
% Y. Tian, T. Zhang, J. Xiao, X. Zhang, and Y. Jin, A coevolutionary
% framework for constrained multi-objective optimization problems, IEEE
% Transactions on Evolutionary Computation, 2020.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    type = Global.ParameterSet(1);

    %% Generate random population
    Population1 = Global.Initialization();
    Population2 = Global.Initialization();
    Fitness1    = CalFitness(Population1.objs,Population1.cons);
    Fitness2    = CalFitness(Population2.objs);

    %% Optimization
    while Global.NotTermination(Population1)
        if type == 1
            MatingPool1 = TournamentSelection(2,Global.N,Fitness1);
            MatingPool2 = TournamentSelection(2,Global.N,Fitness2);
            Offspring1  = GAhalf(Population1(MatingPool1));
            Offspring2  = GAhalf(Population2(MatingPool2));
        elseif type == 2
            MatingPool1 = TournamentSelection(2,2*Global.N,Fitness1);
            MatingPool2 = TournamentSelection(2,2*Global.N,Fitness2);
            Offspring1  = DE(Population1,Population1(MatingPool1(1:end/2)),Population1(MatingPool1(end/2+1:end)));
            Offspring2  = DE(Population2,Population2(MatingPool2(1:end/2)),Population2(MatingPool2(end/2+1:end)));
        end
        [Population1,Fitness1] = EnvironmentalSelection([Population1,Offspring1,Offspring2],Global.N,true);
        [Population2,Fitness2] = EnvironmentalSelection([Population2,Offspring1,Offspring2],Global.N,false);
    end
end