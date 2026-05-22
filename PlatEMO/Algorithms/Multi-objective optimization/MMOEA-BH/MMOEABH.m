classdef MMOEABH < ALGORITHM
% <2025> <multi> <real/integer/label/binary/permutation> <multimodal>
% Multimodal MOEA with block optimization and hybrid clustering
% epsilon --- 0.15 --- The neighborhood radius for DBSCAN
% minpts  ---    5 --- The minimum number of points required to form a dense region

%------------------------------- Reference --------------------------------
% Y. Zhang and W. Hu. Block optimization and switchable hybrid clustering
% for multimodal multiobjective evolutionary optimization with shifted
% local Pareto front. Swarm and Evolutionary Computation, 2025, 98: 102151.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yu Zhang
% If you have any questions, please open an issue in our GitHub repository:
% https://github.com/yuzhang576/MMOEA-BH

    methods
        function main(Algorithm,Problem)
            %% Parameter setting            
            [g, epsilon, minpts] = Algorithm.ParameterSet(10, 0.15, 5);
            
            %% Generate random population
            Population = Problem.Initialization();
            [PBA,LBA]  = deal(cell(1, numel(Population)));                    
            n_PBA      = 5;
            for i = 1:numel(Population)
                [~, temp_PBA,temp_PBA_SCD] = MMOEA_Utils.nd_pccs_sort(Population(i));
                PBA{i}                     = temp_PBA(1:min(numel(temp_PBA),n_PBA));
                PBA_SCD{i}                 = temp_PBA_SCD(1:min(numel(temp_PBA),n_PBA));
            end
            [idx_APC, maxCluster] = MMOEA_Utils.APC(Population.decs);
            for i = 1:maxCluster
                [~, LBA{i}, LBA_SCD{i}] = MMOEA_Utils.nd_pccs_sort(Population(i==idx_APC));
            end
            Archive          = [LBA{1:maxCluster}];
            [gen, count_cso] = deal(1);
            
            %% Optimization
            while Algorithm.NotTerminated(Archive)
                if Problem.FE < Problem.maxFE*0.5 % Stage 1: PSO
                    [Pbest, Nbest] = MMOEA_Utils.rep_selection_pso(Problem, PBA, PBA_SCD, LBA, LBA_SCD, idx_APC, maxCluster);
                    Population     = operator_PSO(Problem, Population, Pbest, Nbest);
                    if mod(gen,g) == 0                        
                        [Archive, LBA, LBA_SCD, PBA, PBA_SCD, idx_APC, maxCluster] = Environmental_Selection('APC', Problem, Archive, Population, LBA, PBA, n_PBA);
                    else                        
                        [Archive, LBA, LBA_SCD, PBA, PBA_SCD, idx_APC, maxCluster] = Environmental_Selection('kmeans', Problem, Archive, Population, LBA, PBA, n_PBA, maxCluster);
                    end
                else                              % Stage 2: CSO
                    if mod(count_cso,2) == 1
                        count_cso          = count_cso + 1;
                        pop_cso            = [];
                        [orig_vol, center] = MMOEA_Utils.get_orig_vol(Problem, maxCluster, LBA);
                        for i = 1 : maxCluster
                            pop                   = LBA{i};
                            [new_up_s, new_low_s] = MMOEA_Utils.get_upper_lower(Problem, pop, orig_vol);
                            [new_up_s, new_low_s] = MMOEA_Utils.clip_upper_lower(new_up_s, new_low_s, center);
                            
                            N                     = ceil(Problem.N/maxCluster);
                            FEs                   = N*25;
                            result                = operator_CSO(Problem, pop, new_low_s, new_up_s, N, FEs);
                            pop_cso               = [pop_cso,result];
                        end
                        if numel(pop_cso) < Problem.N
                            repmat_num = Problem.N - numel(pop_cso);
                            pop_cso    = [pop_cso,pop_cso(randi(numel(pop_cso),repmat_num,1))];
                        end
                        [Archive, LBA, LBA_SCD, PBA, PBA_SCD, idx_APC, maxCluster] = Environmental_Selection('DBSCAN', Problem, Archive, pop_cso(1:Problem.N), LBA, PBA, n_PBA, Problem.N, maxCluster, epsilon, minpts);
                    else
                        count_cso      = count_cso + 1;
                        [Pbest, Nbest] = MMOEA_Utils.rep_selection_pso(Problem, PBA, PBA_SCD, LBA, LBA_SCD, idx_APC, maxCluster);
                        Population     = operator_PSO(Problem, Population, Pbest, Nbest);
                        [Archive, LBA, LBA_SCD, PBA, PBA_SCD, idx_APC, maxCluster] = Environmental_Selection('APC', Problem, Archive, Population, LBA, PBA, n_PBA);
                    end
                end
                gen = gen + 1;
            end
        end
    end
end