classdef BiCo < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained>
% Bidirectional coevolution constrained multiobjective evolutionary algorithm

%------------------------------- Reference --------------------------------
% Z. Liu, B. Wang, and K. Tang, Handling constrained multiobjective
% optimization problems via bidirectional coevolution, IEEE Transactions on
% Cybernetics, 2022, 52(10): 10163-10176.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
       function main(Algorithm,Problem)
           %% Generate random population
           Population = Problem.Initialization();
           ArcPop     = [];
           
           %% Optimization
           while Algorithm.NotTerminated(Population)
               AllPop     = [Population,ArcPop];
               MatingPool = MatingSelection(Population,ArcPop,Problem.N);
               Offspring  = OperatorGA(Problem,MatingPool(1:Problem.N));
               ArcPop     = UpdateArc([AllPop,Offspring],Problem.N);
               Population = EnvironmentalSelection([Population,Offspring],Problem.N);
           end
       end
    end
end