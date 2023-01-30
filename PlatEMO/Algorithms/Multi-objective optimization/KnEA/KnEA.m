classdef KnEA < ALGORITHM
% <many> <real/integer/label/binary/permutation> <constrained/none>
% Knee point driven evolutionary algorithm
% rate --- 0.5 --- Rate of knee points in the population

%------------------------------- Reference --------------------------------
% X. Zhang, Y. Tian, and Y. Jin, A knee point-driven evolutionary algorithm
% for many-objective optimization, IEEE Transactions on Evolutionary
% Computation, 2015, 19(6): 761-776.
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
            rate = Algorithm.ParameterSet(0.5);

            %% Generate random population
            Population = Problem.Initialization();
            FrontNo    = NDSort(Population.objs,Population.cons,inf);
            KneePoints = zeros(1,Problem.N);     % Set of knee points
            r          = -ones(1,2*Problem.N);	% Ratio of size of neighorhood
            t          = -ones(1,2*Problem.N);	% Ratio of knee points

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = MatingSelection(Population.objs,FrontNo,KneePoints);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                Population = [Population,Offspring];
                [FrontNo,MaxFNo]                = NDSort(Population.objs,Population.cons,Problem.N);
                [KneePoints,Distance,r,t]       = FindKneePoints(Population.objs,FrontNo,MaxFNo,r,t,rate);
                [Population,FrontNo,KneePoints] = EnvironmentalSelection(Population,FrontNo,MaxFNo,KneePoints,Distance,Problem.N);      
            end
        end
    end
end