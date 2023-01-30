classdef ICMA < ALGORITHM
% <multi> <real/integer> <constrained>
% Indicator-based constrained multi-objective algorithm

%------------------------------- Reference --------------------------------
% J. Yuan, H. Liu, Y. Ong, and Z. He, Indicator-based evolutionary
% algorithm for solving constrained multi-objective optimization problems,
% IEEE Transactions on Evolutionary Computation, 2022, 26(2): 379-391.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Jiawei Yuan

    methods
        function main(Algorithm,Problem)
            %% Generate the random population
            Population = Problem.Initialization();
            Zmin       = min(Population.objs,[],1);
            Fmin       = min(Population(all(Population.cons<=0,2)).objs,[],1);
            Archive    = Population;
            W          = UniformPoint(Problem.N,Problem.M);
            Ra         = 1;
            
            %% Optimization
            while Algorithm.NotTerminated(Archive)
                Nt = floor(Ra*Problem.N);
                MatingPool = [Population(randsample(Problem.N,Nt)),Archive(randsample(Problem.N,Problem.N-Nt))];
                
                [Mate1,Mate2,Mate3] = Neighbor_Pairing_Strategy(MatingPool,Zmin);
                if rand > 0.5
                    Offspring = OperatorDE(Problem,Mate1,Mate2,Mate3);
                else
                    Offspring = OperatorDE(Problem,Mate1,Mate2,Mate3,{0.5,0.5,0.5,0.75});
                end
                
                Fmin = min([Fmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);
                Zmin = min([Zmin;Offspring.objs],[],1);
                [Population,Archive] = ICMA_Update([Population,Offspring,Archive],Problem.N,W,Zmin,Fmin);
                
                Ra = 1 - Problem.FE/Problem.maxFE;
            end
        end
    end
end