function [CrowdDis, p] = SurvivalScore(front)
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

    distances = zeros(m,m);
    for i=1:m-1
        for j=i+1:m
            % mid point
            mid_point1 = 0.5 * projection(i,:) + projection(j,:) * 0.5;
            t = 1 / sum(mid_point1.^p).^(1/p);
            projected_midpoint1 = mid_point1 * t;

            % Geodesic distance
            distances(i,j) = sum((projection(i,:) - projected_midpoint1).^2)^0.5 + ...
                sum((projection(j,:) - projected_midpoint1).^2)^0.5;

            distances(j, i) = distances(i, j);
        end
    end

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