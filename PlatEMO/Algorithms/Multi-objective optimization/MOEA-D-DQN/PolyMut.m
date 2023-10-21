function Population = PolyMut(Population, lb, ub)
% Polynomial mutation

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    Population = min(max(Population, lb), ub);
    isMut      = rand(size(Population)) < 1 / size(Population, 2);
    u          = rand(size(Population));
    delta1     = (Population - lb) ./ (ub - lb);
    delta2     = (ub - Population) ./ (ub - lb);
    delta      = ((2.*u+(1-2.*u).*(1-delta1).^21).^(1/21)-1).*(u<=0.5) + (1-(2.*(1-u)+(2.*u-1).*(1-delta2).^21).^(1/21)).*(u>0.5);
    Population = Population + isMut .* delta .* (ub - lb);
    Population = min(max(Population, lb), ub);
end