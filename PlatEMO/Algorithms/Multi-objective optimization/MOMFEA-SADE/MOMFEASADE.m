classdef MOMFEASADE < ALGORITHM
% <2022> <multi> <real/integer/label/binary/permutation> <constrained/none> <multitask>
% Multi-objective multifactorial evolutionary algorithm with subspace alignment and adaptive differential evolution

%------------------------------- Reference --------------------------------
% Z. Liang, H. Dong, C. Liu, W. Liang, and Z. Zhu. Evolutionary
% multitasking for multiobjective optimization with subspace alignment and
% adaptive differential evolution. IEEE Transactions on Cybernetics, 2022,
% 52(4): 2096-2109.
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
            RMP   = 1;
            LP    = 30;
            F1    = 0.6;
            F2    = 0.5;
            LCR   = 0.3;
            UCR   = 0.9;
            ProbT = size(Problem.SubM,2);
            for t = 1 : ProbT
                Dec = [rand(Problem.N/2, max(Problem.SubD)),t*ones(Problem.N/2,1)];
                SubPopulation{t} = Problem.Evaluation(Dec,[false*ones(Problem.N/2,1),0*ones(Problem.N/2,1)]);
            end
            STNum  = 3;
            R_Used = []; R_Succ = [];
            ProbT  = size(Problem.SubM,2);
            ProbN  = Problem.N/2;          

            %% Optimization
            Gen = 0;
            while Algorithm.NotTerminated([SubPopulation{:}])
                Gen = Gen+1;
                if Gen <= LP
                    DE_Pro = ones(1, STNum);
                else
                    DE_Pro = sum(R_Succ(Gen - LP:Gen - 1, :) + 1, 1) ./ sum(R_Used(Gen - LP:Gen - 1, :) + 1, 1);
                end                
                Best = getBest(SubPopulation);
                for t = 1 : ProbT
                    % DE Strategies
                    DE_Pool{t} = getDEPool(DE_Pro, length(SubPopulation{t}));
                    % Generation
                    for i = 1 : length(SubPopulation{t})
                        SubPopulation{t}(i).add(1) = false;
                    end
                    offspring = Generation(Problem,SubPopulation, Best, DE_Pool{t}, t,RMP,LP,F1,F2,LCR,UCR);
                    % Selection
                    SubPopulation{t} = [SubPopulation{t}, offspring];
                    rank = NSGA2Sort(SubPopulation{t});
                    SubPopulation{t} = SubPopulation{t}(rank(1:ProbN));
                end
                % DE Strategies Probabilities Updation
                R_Used(Gen, :) = hist([DE_Pool{:}], 1:STNum);
                pop_all = [SubPopulation{:}];
                Adds1   = [];
                Adds2   = [];
                for i = 1 : length(pop_all)
                    Adds1 = [Adds1,pop_all(i).add(1)];                    
                end
                child_idx  = Adds1 == true;
                DE_Succpop = pop_all(child_idx);
                for i = 1 : length(DE_Succpop)                    
                    Adds2 = [Adds2,DE_Succpop(i).add(2)];
                end
                DE_Succ = Adds2;
                R_Succ(Gen, :) = hist(DE_Succ, 1:STNum);
            end
        end
    end
end