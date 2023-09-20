classdef CoMMEA < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <multimodal>
% Coevolutionary multimodal multi-objective evolutionary algorithm
% eps --- 0.2 --- Parameter for quality of the local Pareto front (suggested to be 0 if the problem has no local Pareto front).

%------------------------------- Reference --------------------------------
% W. Li, X. Yao, K. Li, R. Wang, T. Zhang, and  L. Wang, Coevolutionary
% framework for generalized multimodal multi-objective optimization,
% IEEE/CAA Journal of Automatica Sinica, 2023, 10(7): 1544-1556.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Wenhua Li
% p.s. CoMMEA require more function evaluations to receive better
% performance for multimodal multiobjective optimization problems.
    
    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            eps = Algorithm.ParameterSet(0.2);
            
           %% Generate random population
            Population1 = Problem.Initialization();
            Population2 = Problem.Initialization();
            Fitness1    = CalFitness(Population1,'Normal');
            Fitness2    = CalFitness(Population2,'LocalC');
            
           %% Optimization
            while Algorithm.NotTerminated(Population2)
                MatingPool1 = TournamentSelection(2,round(Problem.N/4),Fitness1);
                MatingPool2 = TournamentSelection(2,Problem.N,Fitness2);
                Offspring1  = OperatorGA(Problem,Population1(MatingPool1));
                Offspring2  = OperatorGA(Problem,Population2(MatingPool2));
                
               %% Enviornmental selection
                EvoState = Problem.FE / Problem.maxFE;
                CurEps = max(-log2(1.5*EvoState),eps); % Compute the current epsilon value
                [Population1,Fitness1] = EnvironmentalSelection1([Population1,Offspring1,Offspring2],Problem.N,'Normal',EvoState);
                [Population2,Fitness2] = EnvironmentalSelection2([Population2,Offspring1,Offspring2],Problem.N,'LocalC',CurEps);
            end
        end
    end
end