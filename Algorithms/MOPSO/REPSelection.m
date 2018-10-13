function REP = REPSelection(PopObj,N,div)
% Select one of the particles in REP as the global best position for each
% particle

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    NoP = size(PopObj,1);
    
    %% Calculate the grid location of each solution
    fmax = max(PopObj,[],1);
    fmin = min(PopObj,[],1);
    d    = (fmax-fmin)/div;
    fmin = repmat(fmin,NoP,1);
    d    = repmat(d,NoP,1);
    GLoc = floor((PopObj-fmin)./d);
    GLoc(GLoc>=div) = div - 1;
    GLoc(isnan(GLoc)) = 0;
    
    %% Detect the grid of each solution belongs to
    [~,~,Site] = unique(GLoc,'rows');

    %% Calculate the crowd degree of each grid
    CrowdG = hist(Site,1:max(Site));
    
    %% Calculate the cumulative sum fitnesses of the grids
    Fitness = cumsum(10./CrowdG);
    Fitness = Fitness/max(Fitness);
    
    %% Roulette-wheel selection
    REP = zeros(1,N);
    for i = 1 : length(REP)
        TheGrid = find(rand<=Fitness,1);
        InGrid  = find(Site==TheGrid);
        Temp    = randi(length(InGrid));
        REP(i)  = InGrid(Temp);
    end
end