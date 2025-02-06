function tFrontNo = tNCDSort(PopObj,PopCon,W,fr)
% theta-non-constrained-dominated sorting

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming (email: 20151000334@cug.edu.cn)

    N  = size(PopObj,1);
    NW = size(W,1);
    
    [z,znad] = deal(min(PopObj),max(PopObj));
    [z_c,znad_c] = deal(min(PopCon),max(PopCon));
    %% Normalization
    [PopObj,PopCon] = Normalization(PopObj,PopCon,z,znad,z_c,znad_c);

    %% Calculate the d1 d2 and d3 values for each solution to each weight
    normP  = sqrt(sum(PopObj.^2,2));
    Cosine = 1 - pdist2(PopObj,W,'cosine');
    d1     = repmat(normP,1,size(W,1)).*Cosine;
    d2     = repmat(normP,1,size(W,1)).*sqrt(1-Cosine.^2);
    d3     = sum(max(0,PopCon),2);
    
    %% Clustering
    [~,class] = min(d2,[],2);
    %% Sort
    theta = zeros(1,NW) + 5;
    theta(sum(W>1e-4,2)==1) = 1e6;
    tFrontNo = zeros(1,N);
    for i = 1 : NW
        C = find(class==i);
        [~,rank] = sort(d1(C,i)+theta(i)*d2(C,i)+(1-fr)*d3(i));
        tFrontNo(C(rank)) = 1 : length(C);
    end
end