function Population = PopulationReselection(Grid, GridIndex, N)
    [FrontNo, ~] = NDSort(roundn(Grid.objs, -10), N);
    [~, M] = size(Grid.objs);
    D = pdist2(GridIndex, GridIndex, 'chebychev');
    D(D>1) = 0;
    crowding = sum(D, 2).^(1/M);
    fitness = crowding / max(max(crowding), 1e-10) * (N-1) + 1;
    fitness(FrontNo > 1) = fitness(FrontNo > 1) + max(fitness);
    Population = Grid(RouletteWheelSelection(N, fitness));   
end