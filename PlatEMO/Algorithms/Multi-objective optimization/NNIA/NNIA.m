classdef NNIA < ALGORITHM
% <multi> <real/integer/label/binary/permutation>
% Nondominated neighbor immune algorithm
% nA ---  20 --- Size of active population
% nC --- 100 --- Size of clone population

%------------------------------- Reference --------------------------------
% M. Gong, L. Jiao, H. Du, and L. Bo, Multiobjective immune algorithm with
% nondominated neighbor-based selection, Evolutionary Computation, 2008,
% 16(2): 225-255.
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
            [nA,nC] = Algorithm.ParameterSet(20,100);

            %% Generate random population
            B = Problem.Initialization();               % Antibody population
            D = UpdateDominantPopulation(B,Problem.N);	% Dominant population

            %% Optimization
            while Algorithm.NotTerminated(D)
                A  = D(1:min(nA,length(D)));            % Active population
                C  = Cloning(A,nC);                     % Clone population
                C1 = OperatorGAhalf(Problem,[C,A(randi(length(A),1,length(C)))]);
                D  = UpdateDominantPopulation([D,C1],Problem.N);
            end
        end
    end
end