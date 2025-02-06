classdef CMODEFTR < ALGORITHM
% <2023> <multi> <real/integer> <constrained>
% Constrained multiobjective differential evolution based on the fusion of two rankings

%------------------------------- Reference --------------------------------
% Z. Zeng, X. Zhang, and Z. Hong. A constrained multiobjective differential
% evolution algorithm based on the fusion of two rankings. Information
% Sciences, 2023, 647, 119572.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhiqiang Zeng (email: zhiqiang.zeng@outlook.com)

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                Offspring  = DEgenerator2(Population,Problem);
                t          = Problem.FE/Problem.N;
                MaxGen     = Problem.maxFE/Problem.N;
                a          = 0.5*(1-cos((1-t/MaxGen)*pi));
                Population = EnvironmentalSelection([Population,Offspring],Problem.N,a,Problem);
            end
        end
    end
end