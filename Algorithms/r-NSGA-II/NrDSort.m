function [FrontNo,MaxFNo] = NrDSort(PopObj,nSort,Points,W,delta)
% Do non-r-dominated sorting by efficient non-dominated sort (ENS)

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    FrontNo = inf(1,size(PopObj,1));
    for i = 1 : size(Points,1)
        FrontNo = min(FrontNo,nrdsort(PopObj,Points(i,:),W(i,:),delta));
    end
    MaxFNo = find(cumsum(hist(FrontNo,1:max(FrontNo)))>=min(nSort,length(FrontNo)),1);
    FrontNo(FrontNo>MaxFNo) = inf;
end

function FrontNo = nrdsort(PopObj,g,w,delta)
% Sort the population according to one preferred point

    [PopObj,~,Loc] = unique(PopObj,'rows');
    % Calculate the weighted Euclidean distance of each solution
    Dist       = sqrt(((PopObj-repmat(g,size(PopObj,1),1))./repmat(max(PopObj,[],1)-min(PopObj,[],1),size(PopObj,1),1)).^2*(w/sum(w))');
    DistExtent = max(Dist) - min(Dist);
    % Sort the population based on their Dist values, so that a solution
    % cannot r-dominate the solutions having smaller Dist values than it
    [Dist,rank] = sort(Dist);
    PopObj      = PopObj(rank,:);
    % Non-r-dominated sorting by ENS
    [N,M]   = size(PopObj);
    FrontNo = inf(1,N);
    MaxFNo  = 0;
    while any(FrontNo==inf)
        MaxFNo = MaxFNo + 1;
        for i = 1 : N
            if FrontNo(i) == inf
                Dominated = false;
                for j = i-1 : -1 : 1
                    if FrontNo(j) == MaxFNo
                        m = 2;
                        while m <= M && PopObj(i,m) >= PopObj(j,m)
                            m = m + 1;
                        end
                        Dominated = m > M;
                        if ~Dominated
                            Dominated = (Dist(j)-Dist(i))./DistExtent < -delta;
                        end
                        if Dominated
                            break;
                        end
                    end
                end
                if ~Dominated
                    FrontNo(i) = MaxFNo;
                end
            end
        end
    end
    FrontNo(rank) = FrontNo;
    FrontNo = FrontNo(Loc);
end