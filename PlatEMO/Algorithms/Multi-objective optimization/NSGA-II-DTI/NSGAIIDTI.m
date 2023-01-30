classdef NSGAIIDTI < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained/none> <robust>
% NSGA-II of Deb's type I robust version

%------------------------------- Reference --------------------------------
% K. Deb and H. Gupta, Introducing robustness in multi-objective
% optimization, Evolutionary Computation, 2006, 14(4): 463-494.
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
            [PopObjV,PopConV]    = MeanEffective(Problem,Population);
            [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Problem.N,PopObjV,PopConV);
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
                Offspring  = OperatorGA(Problem,Population(MatingPool),{1,10,1,50});
                [OffObjV,OffConV] = MeanEffective(Problem,Offspring);
                [Population,FrontNo,CrowdDis,PopObjV,PopConV] = EnvironmentalSelection([Population,Offspring],Problem.N,[PopObjV;OffObjV],[PopConV;OffConV]);
            end
        end
    end
end