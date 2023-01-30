classdef BCEIBEA < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% Bi-criterion evolution based IBEA
% kappa --- 0.05 --- Fitness scaling factor

%------------------------------- Reference --------------------------------
% M. Li, S. Yang, and X. Liu, Pareto or non-Pareto: Bi-criterion evolution
% in multiobjective optimization, IEEE Transactions on Evolutionary
% Computation, 2016, 20(5): 645-665.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            kappa = Algorithm.ParameterSet(0.05);

            %% Generate random population
            NPC      = Problem.Initialization();
            [PC,nND] = PCSelection(NPC,Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(PC)
                % PC evolving
                NewPC = Exploration(Problem,PC,NPC,nND,Problem.N);
                % NPC selection
                NPC = EnvironmentalSelection([NPC,NewPC],Problem.N,kappa);
                % NPC evolving
                MatingPool = TournamentSelection(2,Problem.N,-CalFitness(NPC.objs,kappa));
                NewNPC     = OperatorGA(Problem,NPC(MatingPool));
                NPC        = EnvironmentalSelection([NPC,NewNPC],Problem.N,kappa);
                % PC selection
                [PC,nND] = PCSelection([PC,NewNPC,NewPC],Problem.N);
            end
        end
    end
end