classdef tDEACPBI < ALGORITHM
% <2023> <multi/many> <real/integer/label/binary/permutation> <constrained/none>
% Theta-dominance based evolutionary algorithm with CPBI

%------------------------------- Reference --------------------------------
% F. Ming, W. Gong, L. Wang, and L. Gao. A constraint-handling technique
% for decomposition-based constrained many-objective evolutionary
% algorithms. IEEE Transactions on Systems, Man, and Cybernetics: Systems,
% 2023, 53(12): 7783-7793.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming (email: 20151000334@cug.edu.cn)

    methods
        function main(Algorithm,Problem)
            %% Generate the reference points and random population
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            Population    = Problem.Initialization();
            [z,znad]      = deal(min(Population.objs),max(Population.objs));
            [z_c,znad_c]  = deal(min(Population.cons),max(Population.cons));

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = randi(Problem.N,1,Problem.N);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                [Population,z,znad,z_c,znad_c] = EnvironmentalSelection([Population,Offspring],W,Problem.N,z,znad,z_c,znad_c);
            end
        end
    end
end