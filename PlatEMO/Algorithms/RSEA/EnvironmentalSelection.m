function [Population,Range] = EnvironmentalSelection(Global,Population,Range,N)
% The environmental selection of RSEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Non-dominated sorting
    [FrontNO,MaxFNO] = NDSort(Population.objs,N);
    Next             = find(FrontNO<=MaxFNO);
    
    %% Environmental selection
    if any(Range(1,:)==Range(2,:))
        Choose = LastSelection(Global,Population(Next).objs,ismember(Next,find(FrontNO<MaxFNO)),N,ceil(sqrt(N)));
    else
        Choose = LastSelection(Global,(Population(Next).objs-repmat(Range(1,:),length(Next),1))./repmat(Range(2,:)-Range(1,:),length(Next),1),ismember(Next,find(FrontNO<MaxFNO)),N,ceil(sqrt(N))); 
    end
    Population = Population(Next(Choose));
	Range(1,:) = min([Range(1,:);Population.objs],[],1);
    Range(2,:) = max(Population(NDSort(Population.objs,1)==1).objs,[],1);
end

function Choose = LastSelection(Global,PopObj,Choose,N,div)
% Select part of the solutions based on the radar grid

    %% Identify the extreme solutions
    [~,Extreme] = min(repmat(sqrt(sum(PopObj.^2,2)),1,size(PopObj,2)).*sqrt(1-(1-pdist2(PopObj,eye(size(PopObj,2)),'cosine')).^2),[],1); %Calculate the extreme points based on PBI
    Choose      = Choose | ismember(1:size(PopObj,1),Extreme);

    %% Calculate the convergence of each solution
	Con = sum(PopObj.^2,2).^0.5;
    Con = Con./max(Con);
    
    %% Calculate the radar grid of each solution
    [Site,RLoc] = RadarGrid(PopObj,div);
    RDis        = pdist2(RLoc,RLoc);
    RDis(logical(eye(length(RDis)))) = inf;
    CrowdG      = zeros(1,max(Site));
    temp        = tabulate(Site(Choose));
    CrowdG(temp(:,1)) = temp(:,2);

    %% Select N solutions one by one
    while sum(Choose) < N
        remainS  = find(~Choose);
        remainG  = unique(Site(remainS));
        bestG    = CrowdG(remainG) == min(CrowdG(remainG));
        current  = remainS(ismember(Site(remainS),remainG(bestG)));
        r        = 1-(Global.evaluated/Global.evaluation)^2;
        fitness  = Global.M.*r.*Con(current) - min(RDis(current,Choose),[],2);
        [~,best] = min(fitness);
        Choose(current(best))       = true;
        CrowdG(Site(current(best))) = CrowdG(Site(current(best))) + 1;
    end
end