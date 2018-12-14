function [G, nGrids, nDivs] = GridStructure(n, m)
%   GridStructure - Generate a grid deviding the first octant by polar
%   coordinates. Each element of grid is given by (m-1) integer numbers. 
%   The structure is similar to that given by polar coordinates.
%
%   [G, nGrids, nDivs] = GridStructure(n, m)
%
%   Input:
%   n - population size
%   m - number of objectives
%
%   Output:
%   G - grid structure
%   nGrids - number of grids (can be larger than provided population size)
%   nDivs - number of divisions in each right angle
%
%   Example:
%   [G, nGrids, nDivs] = GridStructure(225,3)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Roman Denysiuk

% compute number of divisions
nDivs = ceil( power(n, 1/(m-1)) );

% compute number of grids
nGrids = power(nDivs, m-1);

% initialize grids
G = zeros(nGrids, m-1);

% compute grids
D = m-2;
for j = 0:D
    tmp = repmat(1:nDivs, power(nDivs, j), 1);
    G(:,j+1) = repmat( tmp(:), power(nDivs, D-j), 1);
end

end