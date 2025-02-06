classdef EnvEpsilon
% The environmental selection of SPEA2 (epsilon-constraint method)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        name = 'Env2';  % environemtal selection name
    end

    methods
        function Population = do(obj, Population, N, epsilon)
            %% Calculate the fitness of each solution
            Fitness = AdaFitness(Population.objs,Population.cons,epsilon);

            %% Environmental selection
            Next = Fitness < 1;
            if sum(Next) < N
                [~,Rank] = sort(Fitness);
                Next(Rank(1:N)) = true;
            elseif sum(Next) > N
                Del  = obj.Truncation(Population(Next).objs,sum(Next)-N);
                Temp = find(Next);
                Next(Temp(Del)) = false;
            end
            % Population for next generation
            Population = Population(Next);
            Fitness    = Fitness(Next);

            % Sort the population
            [~,rank] = sort(Fitness);
            Population = Population(rank);
        end

        function Del = Truncation(obj, PopObj,K)
            % Select part of the solutions by truncation
        
            %% Truncation
            Distance = pdist2(PopObj,PopObj);
            Distance(logical(eye(length(Distance)))) = inf;
            Del = false(1,size(PopObj,1));
            while sum(Del) < K
                Remain   = find(~Del);
                Temp     = sort(Distance(Remain,Remain),2);
                [~,Rank] = sortrows(Temp);
                Del(Remain(Rank(1))) = true;
            end
        end
    end
end