function tFrontNo = tNDSort(PopObj,W)
% theta-non-dominated sorting

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming

    N  = size(PopObj,1);
    NW = size(W,1);
    
    [z,znad] = deal(min(PopObj),max(PopObj));
    %% Normalization
    PopObj = Normalization(PopObj,z,znad);

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

function [PopObj,z,znad] = Normalization(PopObj,z,znad)
% Normalize the population and update the ideal point and the nadir point

    [N,M] = size(PopObj);

    %% Update the ideal point
    z = min(z,min(PopObj,[],1));
    
    %% Update the nadir point
    % Identify the extreme points
    W = zeros(M) + 1e-6;
    W(logical(eye(M))) = 1;
    ASF = zeros(N,M);
    for i = 1 : M
        ASF(:,i) = max(abs((PopObj-repmat(z,N,1))./(repmat(znad-z,N,1)))./repmat(W(i,:),N,1),[],2);
    end
    [~,extreme] = min(ASF,[],1);
    % Calculate the intercepts
    Hyperplane = (PopObj(extreme,:)-repmat(z,M,1))\ones(M,1);
    a = (1./Hyperplane)' + z;
    if any(isnan(a)) || any(a<=z)
        a = max(PopObj,[],1);
    end
    znad = a;
    
    %% Normalize the population
    PopObj = (PopObj-repmat(z,N,1))./(repmat(znad-z,N,1));
end