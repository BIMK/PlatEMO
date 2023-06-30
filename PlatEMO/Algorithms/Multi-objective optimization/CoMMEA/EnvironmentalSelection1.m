function [Population,Fitness] = EnvironmentalSelection1(Population,N,Operation,EvoState)
% The environmental selection of CoMMEA

    %% SPEA2 selection
    Fitness = CalFitness(Population,Operation);
    Next = Fitness < 1;
    if sum(Next) < N
        [~,Rank] = sortrows([Fitness' -Crowding(Population.decs)]);
        Population = Population(Rank(1:N));
    else
        Population = Population(Next);
        while length(Population)>N
            dist = sort(pdist2(Population.decs,Population.decs));
            CrowdDis = sum(dist(1:3,:));
            [~,index] = min(CrowdDis);
            Population(index) = [];
        end
    end
    if EvoState<0.5
        Fitness=CalFitness(Population,Operation)-Crowding(Population.decs)';
    else
        Fitness=-Crowding(Population.decs)';
    end
end