function [W,B] = updateWeight(obj,objhat,N)
% Update the weight vectors

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Tomoaki Takagi

    fmin = min(obj,[],1);
    fmax = max(obj,[],1);
    obj = unique([normalize(obj,'range');objhat],'rows');

    %% Select the representative objective vectors
    LpNormD = pdist2(obj,obj,'minkowski',0.5);
    Choose  = false(1,size(obj,1));
    % Select the extreme objective vectors
    [~,extreme]     = min(pdist2(obj,eye(size(obj,2)),'cosine'),[],1);
    Choose(extreme) = true;
    % Greedy inclusion distance-based subset slection
    while sum(Choose) < N
        Remain   = find(~Choose);
        [~, rho] = max(min(LpNormD(Remain,Choose),[],2));
        Choose(Remain(rho)) = true;
    end
    obj = obj(Choose,:);

    %% Update the weight vectors
    W = obj.*repmat(fmax-fmin,N,1);
    B = pdist2(W,W);
    [~,B] = sort(B,2);
    B = B(:,1:ceil(N/10));
    W = W./vecnorm(W,1,2);
end 