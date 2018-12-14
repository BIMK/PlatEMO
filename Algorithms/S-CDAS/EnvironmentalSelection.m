function [Population,FrontNo,CrowdDis] = EnvironmentalSelection(Population,N)
% The environmental selection of S-CDAS

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Non-dominated sorting
    % Conventional non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,N);
    % Sort based on self-controlling dominance area of solutions (S-CDAS)
    for i = 1 : MaxFNo
        FrontNo(FrontNo==i) = FrontNo(FrontNo==i) + SCDASSort(Population(FrontNo==i).objs);
    end
    [~,~,FrontNo] = unique(FrontNo);
    MaxFNo        = find(cumsum(hist(FrontNo,1:max(FrontNo)))>=N,1);
    Next          = FrontNo < MaxFNo;
    
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(Population.objs,FrontNo);
    
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
end

function FrontNo = SCDASSort(PopObj)
% Sort the solutions in a same front by S-CDAS

    [N,M] = size(PopObj);

    %% Translate the solutions
    PopObj = PopObj - min(PopObj,[],1);
    
    %% Create the landmark vectors
    % It should be added but not subtracted a tiny constant value for
    % minimization problem, to avoid the extreme solutions being dominated
    L = diag(max(PopObj,[],1)+1e-2);
    
    %% Check whether each solution dominates the others
    rx   = sqrt(sum(PopObj.^2,2));
    wx   = acos(1-pdist2(PopObj,eye(M),'cosine'));
    lx   = pdist2(PopObj,L);
    Rank = ones(1,N);
    for i = 1 : N
        phix    = repmat(asin(rx(i).*sin(wx(i,:))./lx(i,:)),N,1);
        PopObj1 = repmat(rx,1,M).*sin(wx+phix)./sin(phix);
        dominated       = all(repmat(PopObj1(i,:),N,1)<PopObj1,2);
        Rank(dominated) = Rank(dominated) + 1;
    end
    FrontNo = Rank./(max(Rank)+1);
end