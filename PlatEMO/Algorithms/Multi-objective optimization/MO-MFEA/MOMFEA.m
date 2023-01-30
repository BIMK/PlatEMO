classdef MOMFEA < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained/none> <multitask>
% Multi-objective multifactorial evolutionary algorithm
% rmp --- 1 --- Random mating probability

%------------------------------- Reference --------------------------------
% A. Gupta, Y. Ong, L. Feng, and K. C. Tan, Multiobjective multifactorial
% optimization in evolutionary multitasking, IEEE Transactions on
% Cybernetics, 2017, 47(7): 1652-1665.
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
            rmp = Algorithm.ParameterSet(1);
            
            %% Initialize population
            Population    = Problem.Initialization();
            SubPopulation = Divide(Population,length(Problem.SubD));
            
            %% Optimization
            while Algorithm.NotTerminated([SubPopulation{:}])
                [SubPopulation,Rank] = Sort(Problem,SubPopulation);
                Population           = [SubPopulation{:}];
                ParentPool           = Population(TournamentSelection(2,length(Population),[Rank{:}]));
                SubOffspring         = CreateOff(Problem,ParentPool,rmp);
                for i = 1 : length(Problem.SubD)
                    SubPopulation{i} = EnviSelect([SubPopulation{i},SubOffspring{i}],length(SubPopulation{i}));
                end
            end
        end
    end
end