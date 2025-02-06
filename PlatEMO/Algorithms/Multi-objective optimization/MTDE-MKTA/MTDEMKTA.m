classdef MTDEMKTA < ALGORITHM
% <2024> <multi> <real/integer/label/binary/permutation> <constrained/none> <multitask>
% Multitasking differential evolution with multiple knowledge types and transfer adaptation

%------------------------------- Reference --------------------------------
% Y. Li and W. Gong. Multiobjective multitask optimization with multiple
% knowledge types and transfer adaptation. IEEE Transactions on
% Evolutionary Computation, 2024.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            Tau1  = 0.2;
            Tau2  = 0.1;
            ProbT = size(Problem.SubM,2);
            ProbN = Problem.N/2;

            %% Population initialization
            for t = 1:ProbT
                SubPopulation{t} = [];
                for i = 1:ProbN
                    Dec = [rand(1, max(Problem.SubD)),t];
                    F   = 0.2 + rand();
                    CR  = rand();
                    TR  = rand(); % Knowledge transfer rate
                    KP  = rand(); % Knowledge type proportion
                    Slo = Problem.Evaluation(Dec,[F,CR,TR,KP]);
                    SubPopulation{t} = [SubPopulation{t},Slo];
                end
                [SubPopulation{t},Fitness{t}] = spea2EnvironmentalSelection(SubPopulation{t},ProbN);
                decs = SubPopulation{t}.decs;
                model{t}.mean = mean(decs);
                model{t}.std  = std(decs) +1e-100;
            end
            
            %% Optimization
            while Algorithm.NotTerminated([SubPopulation{:}])
                for t = 1 : ProbT
                    [~,rank{t}] = sort(Fitness{t});
                    [~,rank{t}] = sort(rank{t});
                    decs  = SubPopulation{t}.decs;
                    alpha = 0.5;
                    model{t}.mean = alpha * model{t}.mean + (1 - alpha) * mean(decs);
                    model{t}.std  = alpha * model{t}.std + (1 - alpha) * (std(decs)) +1e-100;
                end
                for t = 1 : ProbT
                    offspring{t} = Generation(Problem,SubPopulation, rank, model, t,Tau1,Tau2);
                end
                for t = 1 : ProbT
                    SubPopulation{t} = [SubPopulation{t}, offspring{t}];
                    [SubPopulation{t},Fitness{t}] = spea2EnvironmentalSelection(SubPopulation{t},ProbN);
                end
            end
        end
    end
end