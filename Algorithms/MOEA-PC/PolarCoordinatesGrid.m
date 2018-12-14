function [idx, theta, rho] = PolarCoordinatesGrid(y, r, nDivs)
%   PolarCoordinatesGrid - Calculate index of grid cell and 
%   polar coordinates for a given objective vector
%
%   [idx, theta, rho] = PolarCoordinatesGrid(y, r, nDivs)
%
%   Input:
%   y - individual's objectives
%   r - reference point
%   nDivs - number of divisions in each right angle
%
%   Output:
%   idx - index of grid cell where individual resides
%   theta - polar angles
%   rho - radius

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Roman Denysiuk

%% calculate polar coordinates

% determine dimensionality
m = numel(r);

% calculate radius
rho = norm(y-r);

% calculate angles
theta = zeros(m-1,1);

for i = 0:m-2
    
    u = (y(m-i) - r(m-i))/rho;
    
    for j = 0:i-1
        u = u/cos(theta(j+1));
    end
    
    u = min( max(u, -1), 1); % ensure feasibility
    
    theta(i+1) = asin(u);
    
    theta(i+1) = min( max(theta(i+1), eps), pi/2-eps);
end

%% calculate grid index based on polar coordinates

% calculate grid coordinates (integer values)
G = ceil(2*nDivs*theta/pi);

% calculate index of pop member using grid coordinates
% formula: idx=G(1)*ndiv^0+(G(2)-1)*ndiv^1+(G(3)-1)*ndiv^2 ...

idx = G(1);
for i = 2:numel(G)
    idx = idx + (G(i)-1)*power(nDivs, i-1);
end

end