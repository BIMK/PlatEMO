function [CrowdDis, p] = SurvivalScore(front)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Annibale Panichella

    [m,n] = size(front);
    CrowdDis = zeros(1,m) ;
    selected = false(1,m);

    front = front - min(front);
    [front, Extreme] = front_normalization(front);

    % set the distance for the extreme solutions
    CrowdDis(Extreme) = Inf;
    selected(Extreme) = true;

    % Newton-Raphson method
    p = ComputeGeometry(front, m, n);

    % project points on computed geometry / shape
    nn = vecnorm(front, p, 2);
    projection = zeros(m, n);
    for i = 1:m
        t = 1 / sum(front(i,:).^p).^(1/p);
        projection(i,:) = front(i,:) * t;
    end

    distances = geodesic_distances(projection, p);
    distances = distances ./ repmat(nn, 1, m);

    neighbors = 2;
    remaining = 1:m;
    remaining = remaining(~selected);
    for i=1:m-sum(selected)-1
        maxim = mink(distances(remaining, selected),neighbors,2);
        [d, index] = max(sum(maxim,2));
        best = remaining(index);
        remaining(index) = [];
        selected(best)=true;
        CrowdDis(1,best) = d;
    end
end

function distances = geodesic_distances(front, p)
[m,~] = size(front);
distances = zeros(m,m);
    for i=1:m-1
        for j=i+1:m
            % mid point
            mid_point1 = 0.5 * front(i,:) + front(j,:) * 0.5;
            t = 1 / sum(mid_point1.^p).^(1/p);
            projected_midpoint1 = mid_point1 * t;

            % Geodesic distance
            distances(i,j) = sum((front(i,:) - projected_midpoint1).^2)^0.5 + ...
                sum((front(j,:) - projected_midpoint1).^2)^0.5;

            distances(j, i) = distances(i, j);
        end
    end
end

function [front, Extreme] = front_normalization(front)
    [m,n] = size(front);

    %% Normalization
    % Detect the extreme points
    Extreme = zeros(1,n);
    w       = zeros(n)+1e-6+eye(n);
    for i = 1 : n
        [~,Extreme(i)] = min(max(front./repmat(w(i,:),m,1),[],2));
    end
    % Calculate the intercepts of the hyperplane constructed by the extreme
    % points and the axes
    Hyperplane = front(Extreme,:)\ones(n,1);
    a = 1./Hyperplane;
    if any(isnan(a))
        a = max(front,[],1)';
    end
    % Normalization
    front = front./repmat(a',m,1);

     % shift the ideal point to the origin
    front = front - min(front);
end 