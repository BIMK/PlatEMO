function [Population,z,znad] = EnvironmentalSelection(Population,W,N,z,znad)
% The environmental selection of theta-DEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,N);
    St = find(FrontNo<=MaxFNo);

    %% Normalization
    [PopObj,z,znad] = Normalization(Population(St).objs,z,znad);
    
    %% theta-non-dominated sorting
    tFrontNo = tNDSort(PopObj,W);
    
    %% Selection
    MaxFNo    = find(cumsum(hist(tFrontNo,1:max(tFrontNo)))>=N,1);
    LastFront = find(tFrontNo==MaxFNo);
    LastFront = LastFront(randperm(length(LastFront)));
    tFrontNo(LastFront(1:sum(tFrontNo<=MaxFNo)-N)) = inf;
    Next      = St(tFrontNo<=MaxFNo);
    % Population for next generation
    Population = Population(Next);
end

function tFrontNo = tNDSort(PopObj,W)
% Do theta-non-dominated sorting

    N  = size(PopObj,1);
    NW = size(W,1);

    %% Calculate the d1 and d2 values for each solution to each weight
    normP  = sqrt(sum(PopObj.^2,2));
    Cosine = 1 - pdist2(PopObj,W,'cosine');
    d1     = repmat(normP,1,size(W,1)).*Cosine;
    d2     = repmat(normP,1,size(W,1)).*sqrt(1-Cosine.^2);
    
    %% Clustering
    [~,class] = min(d2,[],2);
    
    %% Sort
    theta = zeros(1,NW) + 5;
    theta(sum(W>1e-4,2)==1) = 1e6;
    tFrontNo = zeros(1,N);
    for i = 1 : NW
        C = find(class==i);
        [~,rank] = sort(d1(C,i)+theta(i)*d2(C,i));
        tFrontNo(C(rank)) = 1 : length(C);
    end
end