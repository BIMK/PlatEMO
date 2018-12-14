function Z = Adaptive(PopObj,Z,N,interval)
% Addition and deletion of reference points

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    M = size(PopObj,2);
    
    %% Addition of reference points
    rho   = Associate(PopObj,Z);
    old_Z = [];
    while any(rho>=2) && ~isequal(old_Z,Z)
        old_Z = Z;
        for i = find(rho>=2)
            p = repmat(Z(i,:),M,1) - interval/M;
            p(logical(eye(M))) = p(logical(eye(M))) + interval;
            Z = [Z;p];
        end
        Z(any(Z<0,2),:) = [];
        [~,index]       = unique(roundn(Z,-4),'rows','stable');
        Z               = Z(index,:);
        rho = Associate(PopObj,Z);
    end
    
    %% Deletion of reference points
    Z(intersect(N+1:size(Z,1),find(~rho)),:) = [];
end

function rho = Associate(PopObj,Z)
% Associate each solution with one reference point

    %% Calculate the distance of each solution to each reference vector
    NormP    = sqrt(sum(PopObj.^2,2));
    Cosine   = 1 - pdist2(PopObj,Z,'cosine');
    Distance = repmat(NormP,1,size(Z,1)).*sqrt(1-Cosine.^2);
    
    %% Associate each solution with its nearest reference point
    [~,pi] = min(Distance',[],1);
    
    %% Calculate the number of associated solutions of each reference point
    rho = hist(pi,1:size(Z,1));
end