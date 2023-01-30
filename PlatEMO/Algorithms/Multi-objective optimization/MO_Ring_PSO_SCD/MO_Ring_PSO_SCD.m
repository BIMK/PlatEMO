classdef MO_Ring_PSO_SCD < ALGORITHM
% <multi> <real/integer> <multimodal>
% Multiobjective PSO using ring topology and special crowding distance

%------------------------------- Reference --------------------------------
% C. Yue, B. Qu, and J. Liang, A multiobjective particle swarm optimizer 
% using ring topology for solving multimodal multiobjective problems,IEEE 
% Transactions on Evolutionary Computation, 2018, 22(5): 805-817.
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
            %% Initialize parameters
            n_PBA = 5;          % Maximum size of PBA
            n_NBA = 3*n_PBA;	% Maximum size of NBA

            %% Generate random population
            mv   = 0.5*(Problem.upper-Problem.lower);
            Vmin = -mv;
            Vmax = mv;
            ParticleDec = Problem.lower+(Problem.upper-Problem.lower).*rand(Problem.N,Problem.D);
            ParticleVel = Vmin+2.*Vmax.*rand(Problem.N,Problem.D);
            Population  = Problem.Evaluation(ParticleDec,ParticleVel);

            %% Initialize personal best archive PBA and Neighborhood best archive NBA
            PBA = cell(1,Problem.N);
            NBA = cell(1,Problem.N);
            for i = 1:Problem.N
                PBA{i} = Population(i);
                NBA{i} = Population(i);
            end

            %% Optimization
            while Algorithm.NotTerminated(Population)
                NBA = UpdateNBA(NBA,n_NBA,PBA);
                Population = Operator(Problem,Population,PBA,NBA);
                PBA = UpdatePBA(Population,PBA,n_PBA);
                if Problem.FE >= Problem.maxFE
                    tempNBA = [];
                    for i = 1:Problem.N
                        tempNBA = [tempNBA,NBA{i}];
                    end
                    [tempNBA,FrontNo,~] = non_domination_scd_sort(tempNBA,Problem.N);
                    Algorithm.NotTerminated(tempNBA(FrontNo==1));
                end
            end
        end
    end
end