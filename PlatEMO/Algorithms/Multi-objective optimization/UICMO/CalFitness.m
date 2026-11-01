function Fitness = CalFitness(PopObj,PopCon,beta)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yubo Wang (email: 1202111230@cug.edu.cn)

    N        = size(PopObj,1);
    CV       = sum(max(0,PopCon),2);
    Distance = pdist2(PopObj,PopObj);
    Feasible = CV == 0;
    WithinBoundary = any(Distance(:,Feasible) < beta,2);
    Feasible = Feasible | (CV > 0 & WithinBoundary);
    Dominate = false(N);
    for i = 1 : N-1
        for j = i+1 : N
            if Feasible(i) && ~Feasible(j)
                Dominate(i,j) = true;
            elseif ~Feasible(i) && Feasible(j)
                Dominate(j,i) = true;
            elseif ~Feasible(i) && ~Feasible(j)
                if CV(i) < CV(j)
                    Dominate(i,j) = true;
                elseif CV(i) > CV(j)
                    Dominate(j,i) = true;
                end
            else
                k = any(PopObj(i,:)<PopObj(j,:)) - any(PopObj(i,:)>PopObj(j,:));
                if k == 1
                    Dominate(i,j) = true;
                elseif k == -1
                    Dominate(j,i) = true;
                end
            end
        end
    end
    S = sum(Dominate,2);
    R = zeros(1,N);
    for i = 1 : N
        R(i) = sum(S(Dominate(:,i)));
    end
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Distance = sort(Distance,2);
    D        = 1./(Distance(:,floor(sqrt(N)))+2);
    Fitness  = R + D';
end