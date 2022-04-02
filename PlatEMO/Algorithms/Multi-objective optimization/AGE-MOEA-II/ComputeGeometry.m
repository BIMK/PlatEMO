function [p] = ComputeGeometry(front, m, n)

d = pdist2(front, zeros(1,n));
Extreme = FindCornerSolutions(front);
d(Extreme) = Inf;
[~, index] = min(d);

point = front(index,:);
x = NewtonRaphsonMethod(point, 0.001);
    if isnan(x) || x<=0
        p = 1;
    else
        p = abs(x);
    end
end
