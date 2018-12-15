function [FrontNo,MaxFNo] = NDSort_SDR(PopObj,nSort)
% Do non-dominated sorting by strengthened dominance relation (SDR)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N      = size(PopObj,1);
    NormP  = sum(PopObj,2);
    cosine = 1 - pdist2(PopObj,PopObj,'cosine');
    cosine(logical(eye(length(cosine)))) = 0;
    Angle  = acos(cosine);

    temp  = sort(unique(min(Angle,[],2)));
    minA  = temp(min(ceil(N/2),end));
    Theta = max(1,(Angle./minA).^1);
    
    dominate = false(N);
    for i = 1 : N-1
        for j = i+1 : N
            if NormP(i)*Theta(i,j) < NormP(j)
                dominate(i,j) = true;
            elseif NormP(j)*Theta(j,i) < NormP(i)
                dominate(j,i) = true;
            end
        end
    end

    FrontNo = inf(1,N);
    MaxFNo  = 0;
    while sum(FrontNo~=inf) < min(nSort,N)
        MaxFNo  = MaxFNo + 1;
        current = ~any(dominate,1) & FrontNo==inf;
        FrontNo(current)    = MaxFNo;
        dominate(current,:) = false;
    end
end