classdef eMOEA < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% Epsilon multi-objective evolutionary algorithm
% epsilon --- 0.06 --- The parameter in grid location calculation

%------------------------------- Reference --------------------------------
% K. Deb, M. Mohan, and S. Mishra, Towards a quick computation of
% well-spread Pareto-optimal solutions, Proceedings of the International
% Conference on Evolutionary Multi-Criterion Optimization, 2003, 222-236.
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
            epsilon = Algorithm.ParameterSet(0.06);

            %% Generate random population
            Population = Problem.Initialization();
            PopGrid    = floor((Population.objs-repmat(min(Population.objs,[],1),Problem.N,1))/epsilon);
            eFrontNO   = NDSort(PopGrid,1);
            Archive    = Population(eFrontNO==1);

            %% Optimization
            while Algorithm.NotTerminated(Archive)
                for i = 1 : Problem.N
                    k    = randperm(Problem.N,2);
                    Domi = any(Population(k(1)).obj<Population(k(2)).obj) - any(Population(k(1)).obj>Population(k(2)).obj);
                    p    = k((Domi==-1)+1);
                    q    = randi(length(Archive));
                    Offspring  = OperatorGAhalf(Problem,[Population(p),Archive(q)]);
                    Population = UpdatePopulation(Population,Offspring);
                    Archive    = UpdateArchive(Archive,Offspring,epsilon);
                end
            end
        end
    end
end