function Fitness = FitnessSingle(Population)
%FitnessSingle - Fitness calculation for single-objective optimization.
%
%   Fit = FitnessSingle(P) calculates the fitness value of each solution in
%   P for single-objective optimization, where both the objective value and
%   constraint violation are considered.
%
%   Example:
%       Fitness = FitnessSingle(Population)

%------------------------------- Reference --------------------------------
% K. Deb, An efficient constraint handling method for genetic algorithms,
% Computer Methods in Applied Mechanics and Engineering, 2000, 186(2-4):
% 311-338.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    PopCon   = sum(max(0,Population.cons),2);
    Feasible = PopCon <= 0;
    Fitness  = Feasible.*Population.objs + ~Feasible.*(PopCon+1e10);
end