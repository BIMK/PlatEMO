classdef AGEA < ALGORITHM
% <multi> <real/binary/permutation>
% Adaptive grid based evolutionary algorithm
% limit --- 1 --- limit output solutions number to N (1. limit  0. no limit)
% div --- 10 --- initial number of grid divisions

%------------------------------- Reference --------------------------------
% Z. Liu, F. Han, Q. Ling, H. Han, J. Jiang, and Q. Liu, A multi-objective
% evolutionary algorithm based on a grid with adaptive divisions for
% multi-objective optimization with irregular Pareto fronts, Applied Soft 
% Computing, 2025, 176: 113106.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
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
            [limit, div]= Algorithm.ParameterSet(1, 10);
            %% Generate random population
            Population = Problem.Initialization(); 
            Grid   = Population;  
            zmin = min(Population.objs,[],1);
            zmax = max(max(Population.objs, [], 1), zmin + 1e-10);
            gmax = zmax;            
            while Algorithm.NotTerminated(Grid)
            %% Optimization
                MatingPool = randperm(length(Population));
                Offspring = OperatorGA(Problem, Population(MatingPool));
                Population = [Grid, Offspring];
                [FrontNo, ~] = NDSort(roundn(Population.objs, -10), Problem.N);              
                Grid = Population(FrontNo ~= inf);
                NDPop = Population(FrontNo == 1);
                zmin = min(zmin, min(Grid.objs,[],1));
                gmax = GridStabilization(NDPop, zmin, gmax, div);               
                [Grid, GridIndex, div] = GridAdaptiveAdjustment(Grid, zmin, gmax, div, Problem.N);              
                Population = PopulationReselection(Grid, GridIndex, Problem.N);                                    
              %% limit output solutions number to N
                if limit ~= 0 && length(Grid) > Problem.N
                    if Problem.FE >= Problem.maxFE
                        Grid = EnvironmentalSelection_SPEA2(Grid, Problem.N);
                    end
                end
            end
        end
    end
end


