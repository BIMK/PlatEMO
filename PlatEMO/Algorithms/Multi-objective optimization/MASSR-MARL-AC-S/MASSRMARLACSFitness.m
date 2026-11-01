function Fitness = MASSRMARLACSFitness(PopObj)
%MASSRMARLACSFitness - Calculate SPEA2 strength and density fitness.
%
%   Raw fitness sums the strengths of dominating solutions. Add the density
%   estimate from the floor(sqrt(N))-th nearest objective-space neighbor.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N = size(PopObj,1);

    % Pairwise Pareto dominance
    Dominate = false(N);
    for i = 1 : N-1
        for j = i+1 : N
            k = any(PopObj(i,:)<PopObj(j,:)) - any(PopObj(i,:)>PopObj(j,:));
            if k == 1
                Dominate(i,j) = true;
            elseif k == -1
                Dominate(j,i) = true;
            end
        end
    end

    % Strength and raw fitness
    S = sum(Dominate,2);

    R = zeros(1,N);
    for i = 1 : N
        R(i) = sum(S(Dominate(:,i)));
    end

    % Nearest-neighbor density
    Distance                                 = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Distance                                 = sort(Distance,2);
    D                                        = 1./(Distance(:,floor(sqrt(N)))+2);

    Fitness = R + D';
end