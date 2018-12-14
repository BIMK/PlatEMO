function RgFrontNO = MaximumRanking(PopObj,W,K)
% Get the front of each solution by calculating the global rank

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N  = size(PopObj,1);
    NW = size(W,1);

    %% Find the K nearest weight vectors of each solution
    [~,rank] = sort(1-pdist2(PopObj,W,'cosine'),2,'descend');
    L        = rank(:,1:K);
    
    %% Calculate the modified Tchebycheff function value
	g = inf(N,NW);
    for i = 1 : N
        g(i,L(i,:)) = max(repmat(PopObj(i,:),K,1)./W(L(i,:),:),[],2)';
    end

    %% Calculate the global rank of each solution
    [~,rank]  = sort(g,1);
    [~,r]     = sort(rank,1);
    g(g~=inf) = 1;
    Rg        = min(r.*g,[],2);
    
    %% Get the front of each solution
    [~,~,RgFrontNO] = unique(Rg);
end