classdef MFEAII < ALGORITHM
% <single> <real/integer/label/binary/permutation> <large/none> <multitask>
% Multifactorial evolutionary algorithm II

%------------------------------- Reference --------------------------------
% K. K. Bali, Y. Ong, A. Gupta, and P. S. Tan, Multifactorial evolutionary
% algorithm with online transfer parameter estimation: MFEA-II, IEEE
% Transactions on Evolutionary Computation, 2020, 24(1): 69-83.
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
                RMP          = learnRMP(Problem,SubPopulation);
                SubParent    = SelectHalf(SubPopulation);
                SubOffspring = Child(Problem,SubParent,RMP);
                for i = 1 : length(Problem.SubD)
                    SubPopulation{i}= [SubPopulation{i},SubOffspring{i}];
                end
                SubPopulation = SelectPop(SubPopulation,Problem);
            end
        end
    end
end