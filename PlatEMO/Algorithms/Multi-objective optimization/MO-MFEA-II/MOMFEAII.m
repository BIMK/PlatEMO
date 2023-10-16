classdef MOMFEAII < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained/none> <multitask>
% Multi-objective multifactorial evolutionary algorithm II

%------------------------------- Reference --------------------------------
% K. K. Bali, A. Gupta, Y. Ong, and P. S. Tan, Cognizant multitasking in
% multiobjective multifactorial evolution: MO-MFEA-II, IEEE Transactions on
% Cybernetics, 2021, 51(4): 1784-1796.
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
            %% Initialize population
            Population    = Problem.Initialization();
            SubPopulation = Divide(Population,length(Problem.SubD));
            
            %% Optimization
            while Algorithm.NotTerminated([SubPopulation{:}])
                RMP                  = learnRMP(Problem,SubPopulation);
                [SubPopulation,Rank] = Sort(Problem,SubPopulation);
                Population           = [SubPopulation{:}];
                ParentPool           = Population(TournamentSelection(2,length(Population),[Rank{:}]));
                SubOffspring         = CreateOff(Problem,ParentPool,SubPopulation,RMP);
                for i = 1 : length(Problem.SubD)
                    SubPopulation{i} = EnviSelect([SubPopulation{i},SubOffspring{i}],length(SubPopulation{i}));
                end
            end
        end
    end
end
