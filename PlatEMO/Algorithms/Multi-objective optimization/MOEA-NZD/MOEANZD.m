classdef MOEANZD < ALGORITHM
% <2024> <multi/many> <real> <large/none> <constrained/none> <sparse>
% Multi-objective evolutionary algorithm with nonzero detection

%------------------------------- Reference --------------------------------
% X. Wang, R. Cheng, and Y. Jin. Sparse large-scale multiobjective
% optimization by identifying nonzero decision variables. IEEE Transactions
% on Systems, Man, and Cybernetics: Systems, 2024, 54(10): 6280-6292.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Xiangyu Wang (email: xiangyu.wang@uni-bielefeld.de)

    methods
        function main(Algorithm,Problem)
            %% Generate the reference points and random population
            [Z,Problem.N] = UniformPoint(Problem.N,Problem.M);
            Mask = zeros(Problem.N, Problem.D); % Generate population with all zeros.
            Population    = Problem.Evaluation(Mask);
            Zmin          = min(Population(all(Population.cons<=0,2)).objs,[],1);

            %% Optimization
            % predefined parameter
            Re0 = false(1, Problem.D);
            Re1 = true(1, Problem.D);
            t   = 1;
            t1  = 1;
            evaluation  = 0;
            Demarcation = floor(Problem.maxFE * 0.7 / Problem.N);
            intervel    = floor(Problem.maxFE / Problem.N / 20);
            while Algorithm.NotTerminated(Population)
                evaluation = evaluation + 1;
                MatingPool = TournamentSelection(2,Problem.N,sum(max(0,Population.cons),2));
                NextStep   = mean(abs(Re1-Re0)) >= (1/Problem.D);
                if  evaluation < Demarcation
                    if NextStep && mod(evaluation,intervel) == 0
                        if t == 1
                            t = t + 1;
                        else
                            Re0 = Re1;
                            t   = t + 1;
                        end
                        [Re1, PopulationM] = DimJud(Population(MatingPool).decs, Problem.upper, Problem.lower, Re0);
                    else
                        PopulationM = Population(MatingPool).decs;
                    end
                    Offspring = OperatorGA(Problem, PopulationM, {1,20,1,1});
                    Offspring = Problem.Evaluation(Offspring);
                else
                    Offspring0 = Population(MatingPool).decs;
                    if t1 == 1 || mod(evaluation,intervel) == 0
                        [Offspring,Non0_index] = DimJud0(Offspring0, Problem);
                        t1 = t1 + 1;
                    else
                        Offspring = GA0(Offspring0,Non0_index, Problem); 
                    end
                    Offspring  = Problem.Evaluation(Offspring);
                    Offspring0 = [];
                end
                Zmin       = min([Zmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);
                Population = EnvironmentalSelection([Population,Offspring],Problem.N,Z,Zmin);
            end
        end
    end
end