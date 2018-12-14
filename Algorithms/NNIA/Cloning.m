function C = Cloning(A,nC)
% Proportional cloning

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(A.objs);
    if all(CrowdDis==inf)
        CrowdDis = ones(size(CrowdDis));
    else
        CrowdDis(CrowdDis==inf) = 2*max(CrowdDis(CrowdDis~=inf));
    end

    %% Clone
    q = ceil(nC*CrowdDis/sum(CrowdDis));
    C = [];
    for i = 1 : length(A)
        C = [C,repmat(A(i),1,q(i))];
    end
end

function CrowdDis = CrowdingDistance(PopObj)
% Calculate the crowding distance of each solution in the same front

    [N,M]    = size(PopObj);
    
    CrowdDis = zeros(1,N);
    Fmax     = max(PopObj,[],1);
    Fmin     = min(PopObj,[],1);
    for i = 1 : M
        [~,rank] = sortrows(PopObj(:,i));
        CrowdDis(rank(1))   = inf;
        CrowdDis(rank(end)) = inf;
        for j = 2 : N-1
            CrowdDis(rank(j)) = CrowdDis(rank(j))+(PopObj(rank(j+1),i)-PopObj(rank(j-1),i))/(Fmax(i)-Fmin(i));
        end
    end
end