classdef DGEA < ALGORITHM
% <multi/many> <real/integer> <large/none>
% Direction guided evolutionary algorithm
% operation ---   1 --- Operation of the environmental selection
% RefNo     ---  10 --- Number of reference vectors for offspring generation

%------------------------------- Reference --------------------------------
% C. He, R. Cheng, and D. Yazdani, Adaptive offspring generation for
% evolutionary large-scale multiobjective optimization, IEEE Transactions
% on System, Man, and Cybernetics: Systems, 2022, 52(2): 786-798.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            [operation, RefNo]= Algorithm.ParameterSet(1, 10);
            [V,Problem.N] = UniformPoint(Problem.N,Problem.M);
            Population    = Problem.Initialization();
            Offspring     = Problem.Initialization();
            Arc           = Population;

            %% Optimization
            while Algorithm.NotTerminated(Arc)
                [Population,FrontNo] = PreSelection([Population,Offspring],V,(Problem.FE/Problem.maxFE)^2,RefNo);
                Offspring = DirectionReproduction(Problem,Population,FrontNo,RefNo);
                switch operation
                    case 1 
                        Arc = subRVEA([Arc,Offspring],V,(Problem.FE/Problem.maxFE)^2);
                    case 2
                        Arc = subNSGAII([Arc,Offspring],Problem.N);
                    case 3 
                        Arc = subIBEA([Arc,Offspring],Problem.N,0.05);
                    case 4 
                        Arc = subSPEA2([Arc,Offspring],Problem.N);
                    otherwise
                        Arc = [Population,Offspring];	% Without selection
                end
            end
        end
    end
end