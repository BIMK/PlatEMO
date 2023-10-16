classdef DSPCMDE < ALGORITHM
% <multi> <real/integer> <constrained>
% Dynamic selection preference-assisted constrained multiobjective differential evolution

%------------------------------- Reference --------------------------------
% K. Yu, J. Liang, B. Qu, Y. Luo, and C. Yue, Dynamic selection
% preference-assisted constrained multiobjective differential evolution,
% IEEE Transactions on Systems, Man, and Cybernetics: Systems, 2022, 52(5):
% 2954-2965.
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
                Offspring  = DEgenerator2(Problem,Population);
                a          = 0.5*(1-cos((1-Problem.FE/Problem.maxFE)*pi));                
                Population = EnvironmentalSelection([Population,Offspring],Problem.N,a);
            end
        end
    end
end