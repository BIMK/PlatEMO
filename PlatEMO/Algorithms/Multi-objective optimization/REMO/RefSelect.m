function Ref = RefSelect(Population,k)
% Reference solutions selection by RSEA strategy

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    k      = min(k,length(Population));
    PopObj = Population.objs;
	[FrontNO,MaxFNO] = NDSort(PopObj,k);
    Next = find(FrontNO<=MaxFNO);
    Pmin = min(PopObj,[],1) + 1e-6;
    Pmax = max(PopObj,[],1);
    if Pmax > Pmin
        PopObj = (PopObj-repmat(Pmin,size(PopObj,1),1))./repmat(Pmax-Pmin,size(PopObj,1),1);
    end
    
    %% Environmental selection
    Choose = LastSelection(PopObj(Next,:),ismember(Next,find(FrontNO<MaxFNO)),ceil(sqrt(k)),k);
    Ref    = Population(Next(Choose));
end
    
function Choose = LastSelection(PopObj,Choose,div,k)
% Select part of the solutions based on the radar grid
    
    %% Identify the extreme solutions
	[~,Extreme] = min(sqrt(sum(PopObj.^2,2)).*sqrt(1-(1-pdist2(PopObj,ones(1,size(PopObj,2)),'cosine')).^2),[],1); %Calculate the extreme points based on PBI
    Choose      = Choose | ismember(1:size(PopObj,1),Extreme);

    %% Calculate the convergence of each solution
	Con = sum(PopObj.^1,2).^1;
    Con = Con./max(Con);
    
    %% Calculate the radar grid of each solution
    [Site,RLoc] = RadarGrid(PopObj,div);
    RDis        = pdist2(RLoc,RLoc);
    RDis(logical(eye(length(RDis)))) = inf;
    CrowdG      = zeros(1,max(Site));
    temp        = tabulate(Site(Choose));
    CrowdG(temp(:,1)) = temp(:,2);

    %% Select k solutions
    while sum(Choose) < k
        % Delete outline solutions
        remainS  = find(~Choose);
        remainG  = unique(Site(remainS));
        bestG    = CrowdG(remainG) == min(CrowdG(remainG));
        current  = remainS(ismember(Site(remainS),remainG(bestG)));
        fitness  = 0.1.*size(PopObj,2).*Con(current) - min(RDis(current,Choose),[],2); % - 0.1.* min(Dis(current,Choose),[],2);
        [~,best] = min(fitness);
        Choose(current(best))       = true;
        CrowdG(Site(current(best))) = CrowdG(Site(current(best))) + 1;
    end
end   

function [Site,RLoc] = RadarGrid(P,div)

	[N,M] = size(P);
     
    %% Calculate the radar coordinate of each solution
    theta     = 0 : 2*pi/M : 2*pi/M*(M-1);
    RLoc(:,1) = sum(P.*repmat(cos(theta),N,1),2)./sum(P,2);
    RLoc(:,2) = sum(P.*repmat(sin(theta),N,1),2)./sum(P,2);
    RLoc      = (RLoc+1)/2;
    YL        = min(RLoc,[],1);                             % Lower bounary of the transferred points
    YU        = max(RLoc,[],1);                             % Upper bounary of the transferred points  
    NRLoc     = (RLoc-repmat(YL,N,1))./repmat(YU-YL,N,1);	% Normalized points
    
    %% Identify the index of grid of each solution
    GLoc            = floor(NRLoc.*div);
    GLoc(GLoc>=div) = div - 1;
    UniqueGLoc      = sortrows(unique(GLoc,'rows'));
    [~,Site]        = ismember(GLoc,UniqueGLoc,'rows');
end