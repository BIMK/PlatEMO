classdef HEA < ALGORITHM
% <multi/many> <real/binary/permutation>
% Hyper-dominance based evolutionary algorithm
% MaxT --- 0.05 --- The maximum value of tolerance

%------------------------------- Reference --------------------------------
% Z. Liu, F. Han, Q. Ling, H. Han, and J. Jiang, A many-objective 
% optimization evolutionary algorithm based on hyper-dominance degree, 
% Swarm and Evolutionary Computation, 2023, 83: 101411.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhe Liu

    methods
        function main(Algorithm,Problem)
            
            %% Parameter setting
            MaxT = Algorithm.ParameterSet(0.05);
            [W, Problem.N] = UniformPoint(Problem.N,Problem.M);
            T = 0;
            step = MaxT/(Problem.maxFE / Problem.N);
            
            %% Generate random population
            Population          = Problem.Initialization(); 
            Solutions =  Population;
            zmin = min(Population.objs,[],1);
            zmax = max(max(Population.objs, [], 1), zmin + 1e-6);
            
            %% Optimization
            while Algorithm.NotTerminated(Solutions)       
                MatingPool = randperm(length(Population));
                Offsprings = OperatorGA(Problem, Population(MatingPool));
                Population = [Offsprings, Solutions];  
                zmin = min(zmin, min(Population.objs,[],1));          
                [non_dominated, hd] = DominationCal_HEA(Population.objs, zmin, zmax, T);               
                Population = Population(non_dominated == 1);
                hd = hd(non_dominated == 1);
                zmax = max(max(Population.objs, [], 1), zmin + 1e-6);
                T = T + step;
                [Solutions, hd] = EnvironmentalSelection_HEA(Population, hd, zmin, zmax, Problem.N, W); 
                
                %% Population reselection strategy
                Population = Solutions;
                r = unidrnd(Problem.N,[1,Problem.N]);
                for i = 1: Problem.N
                    if hd(i) < hd(r(i))
                        Population(i) = Solutions(r(i));
                    end
                end
            end
        end
    end
end