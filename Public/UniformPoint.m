function [W,N] = UniformPoint(N,M)
%UniformPoint - Generate a set of uniformly distributed points on the unit
%hyperplane.
%
%   [W,N] = UniformPoint(N,M) returns approximately N uniformly distributed
%   points with M objectives on the unit hyperplane.
%
%   Due to the requirement of uniform distribution, the number of points
%   cannot be arbitrary, and the number of points in W may be slightly
%   smaller than the predefined size N.
%
%   Example:
%       [W,N] = UniformPoint(275,10)

%------------------------------- Reference --------------------------------
% [1] I. Das and J. E. Dennis, Normal-boundary intersection: A new method
% for generating the Pareto surface in nonlinear multicriteria optimization
% problems, SIAM Journal on Optimization, 1998, 8(3): 631-657.
% [2] K. Deb and H. Jain, An evolutionary many-objective optimization
% algorithm using reference-point based non-dominated sorting approach,
% part I: Solving problems with box constraints, IEEE Transactions on
% Evolutionary Computation, 2014, 18(4): 577-601.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    H1 = 1;
    while nchoosek(H1+M,M-1) <= N
        H1 = H1 + 1;
    end
    W = nchoosek(1:H1+M-1,M-1) - repmat(0:M-2,nchoosek(H1+M-1,M-1),1) - 1;
    W = ([W,zeros(size(W,1),1)+H1]-[zeros(size(W,1),1),W])/H1;
    if H1 < M
        H2 = 0;
        while nchoosek(H1+M-1,M-1)+nchoosek(H2+M,M-1) <= N
            H2 = H2 + 1;
        end
        if H2 > 0
            W2 = nchoosek(1:H2+M-1,M-1) - repmat(0:M-2,nchoosek(H2+M-1,M-1),1) - 1;
            W2 = ([W2,zeros(size(W2,1),1)+H2]-[zeros(size(W2,1),1),W2])/H2;
            W  = [W;W2/2+1/(2*M)];
        end
    end
    W = max(W,1e-6);
    N = size(W,1);
end