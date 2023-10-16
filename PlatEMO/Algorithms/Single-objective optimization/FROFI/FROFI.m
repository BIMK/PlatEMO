classdef FROFI < ALGORITHM
% <single> <real/integer> <large/none> <constrained>
% Feasibility rule with the incorporation of objective function information

%------------------------------- Reference --------------------------------
% Y. Wang, B. Wang, H. Li, G. G. Yen, Incorporating objective function
% information into the feasibility rule for constrained evolutionary
% optimization, IEEE Transactions on Cybernetics, 2015, 46(12): 1-15.
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
            %% Generate random population
            Population = Problem.Initialization();
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                Offspring  = Operator(Problem,Population);
                Population = EnvironmentalSelection(Population,Offspring);
                Population = Mutation(Problem,Population);
            end
        end
    end
end