function [p] = ComputeGeometry(front, m, n)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Annibale Panichella

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