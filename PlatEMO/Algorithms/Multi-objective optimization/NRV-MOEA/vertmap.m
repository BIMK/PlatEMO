function [mapPop,Hyperplane] = vertmap(ArcObj,PopObj,Hyperplane_bp,Global)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [N,M]  = size(PopObj);
    mapPop = zeros(N,M);
    % Find the extreme points
    [~,Rank]   = sort(ArcObj,'descend');
    Extreme    = zeros(1,M);
    Extreme(1) = Rank(1,1);
    for j = 2 : length(Extreme)
        k = 1;
        Extreme(j) = Rank(k,j);
        while ismember(Extreme(j),Extreme(1:j-1)) && k < size(Rank,1)
            k = k+1;
            Extreme(j) = Rank(k,j);
        end
    end
    % Calculate the hyperplane
    try
        if size(ArcObj,1)>=M
            Hyperplane = ArcObj(Extreme,:)\ones(length(Extreme),1);
        else
            Hyperplane = PopObj(Extreme,:)\ones(length(Extreme),1);
        end
    catch
        Hyperplane = ones(M,1);
    end
    % Calculate the map point of each solution to the hyperplane
    if sum(isinf(Hyperplane)) > 0
        Hyperplane = ones(M,1);
    elseif sum(isnan(Hyperplane)) > 0
        Hyperplane = ones(M,1);
    end
    for i = 1 : N
        p  = PopObj(i,:);
        t1 = sum(p'.*Hyperplane)-1;
        t2 = sum(Hyperplane.^2);
        for m = 1 : M
            mapPop(i,m) = (-Hyperplane(m)*(t1-Hyperplane(m)*p(m))+p(m)*(t2-Hyperplane(m)^2))/t2;
        end
    end
end