function Population = EnvironmentalSelection(Population,N,div)
% The environmental selection of PESA-II

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Select the non-dominated solutions
    Next = NDSort(Population.objs,1) == 1;

    %% Delete the solutions to make a total of N solutions be chosen
    if sum(Next) > N
        Del  = Delete(Population(Next).objs,sum(Next)-N,div);
        Temp = find(Next);
        Next(Temp(Del)) = false;
    end
    % Population for next generation
    Population = Population(Next);
end

function Del = Delete(PopObj,K,div)
% Delete part of the non-dominated solutions by the grid

    N = size(PopObj,1);

    %% Calculate the grid location of each solution
    fmax = max(PopObj,[],1);
    fmin = min(PopObj,[],1);
    d    = (fmax-fmin)/div;
    GLoc = floor((PopObj-repmat(fmin,N,1))./repmat(d,N,1));
    GLoc(GLoc>=div)   = div - 1;
    GLoc(isnan(GLoc)) = 0;

    %% Calculate the crowding degree of each grid
    [~,~,Site] = unique(GLoc,'rows');
    CrowdG     = hist(Site,1:max(Site));

    %% Delete K solutions
    Del = false(1,N);
    while sum(Del) < K
        % Select the most crowded grid
        maxGrid = find(CrowdG==max(CrowdG));
        Temp    = randi(length(maxGrid));
        Grid    = maxGrid(Temp);
        % And delete one solution randomly from the grid
        InGrid  = find(Site==Grid);
        Temp    = randi(length(InGrid));
        p       = InGrid(Temp);
        Del(p)  = true;
        Site(p) = NaN;
        CrowdG(Grid) = CrowdG(Grid) - 1;
    end
end