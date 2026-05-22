function [Population, FrontNo, CrowdDis] = EnvironmentalSelectionGW(Vectors, Population,Utop, Nadir, nsort, ro, eps)
    %% Non-dominated sorting
    [FrontNo, MaxFNo] = GWASFGASort(Vectors, Population.objs,Utop, Nadir, nsort, ro, eps);
    Next = FrontNo < MaxFNo;

    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(Population.objs, FrontNo);

    %% Select the solutions in the last front by their crowding distances
    Last = find(FrontNo == MaxFNo);
    [~, Rank] = sort(CrowdDis(Last), 'descend');
    numSelected = min(nsort - sum(Next), numel(Last));  % Avoid selecting more than available
    Next(Last(Rank(1:numSelected))) = true;

    %% Population for next generation
    Population = Population(Next);
    FrontNo = FrontNo(Next);
    CrowdDis = CrowdDis(Next);
end
