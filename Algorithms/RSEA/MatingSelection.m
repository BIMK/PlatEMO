function MatingPool = MatingSelection(PopObj,Range,N)
% The mating selection of RSEA

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

	%% Calculate the convergence of each solution
    PopObj = (PopObj-repmat(Range(1,:),size(PopObj,1),1))./repmat(Range(2,:)-Range(1,:),size(PopObj,1),1);
    Con = sum(PopObj.^2,2).^0.5;

    %% Calculate the radar grid of each solution
    [Site,~] = RadarGrid(PopObj,ceil(sqrt(size(PopObj,1)))); 
    temp     = tabulate(Site);
    CrowdG   = temp(:,2);
    
    %% Binary tournament selection
    MatingPool = zeros(1,ceil(N/2)*2);
    grids      = TournamentSelection(2,length(MatingPool),CrowdG);
    for i = 1 : length(MatingPool)
        current = find(Site==grids(i));
        if isempty(current)
             MatingPool(i) = randi(size(PopObj,1),1);
        else
            parents       = current(randi(length(current),1,2));
            [~,best]      = min(Con(parents));
            MatingPool(i) = parents(best);
        end
    end
end