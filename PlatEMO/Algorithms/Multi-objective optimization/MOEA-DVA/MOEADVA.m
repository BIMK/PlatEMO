classdef MOEADVA < ALGORITHM
% <multi> <real/integer> <large>
% Multi-objective evolutionary algorithm based on decision variable
% analyses
% NCA --- 20 --- The number of sampling solutions in control variable analysis
% NIA ---  6 --- The maximum number of tries required to judge the interaction

%------------------------------- Reference --------------------------------
% X. Ma, F. Liu, Y. Qi, X. Wang, L. Li, L. Jiao, M. Yin, and M. Gong, A
% multiobjective evolutionary algorithm based on decision variable analyses
% for multiobjective optimization problems with large-scale variables, IEEE
% Transactions Evolutionary Computation, 2016, 20(2): 275-298.
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
            [NCA,NIA] = Algorithm.ParameterSet(20,6);

            %% Control variable analysis
            [DiverIndexes,ConverIndexes] = ControlVariableAnalysis(Problem,NCA);

            %% Dividing distance variables based on two variable analyses
            [Subcomponents,Population] = DividingDistanceVariables(Problem,NIA,DiverIndexes,ConverIndexes);

            %% Calculate the neighbours of each individual
            PopDec = Population.decs;
            Dis    = pdist2(PopDec(:,DiverIndexes),PopDec(:,DiverIndexes));
            Dis(logical(eye(length(Dis)))) = inf;
            [~,Neighbour] = sort(Dis,2);
            Neighbour     = Neighbour(:,1:ceil(Problem.N/10));

            %% Subcomponent optimization
            if Problem.M == 2; threshold = 0.01; else threshold = 0.03; end
            while Algorithm.NotTerminated(Population)
                for i = 1 : length(Subcomponents)
                    drawnow('limitrate');
                    Population = SubcomponentOptimizer(Problem,Population,Neighbour,Subcomponents{i});
                end
            end
        end
    end
end