classdef GWO < ALGORITHM
% <2014> <single> <real/integer> <large/none> <constrained/none>
% Grey wolf optimizer

%------------------------------- Reference --------------------------------
% S. Mirjalili, S. M. Mirjalili, and A. Lewis. Grey wolf optimizer.
% Advances in Engineering Software, 2014, 69: 46-61.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                [~,r]  = sort(FitnessSingle(Population));
                alpha  = repmat(Population(r(1)).dec,Problem.N,1);
                beta   = repmat(Population(r(2)).dec,Problem.N,1);
                delta  = repmat(Population(r(3)).dec,Problem.N,1);
                Dalpha = abs(2*rand(Problem.N,Problem.D).*alpha-Population.decs);
                Dbeta  = abs(2*rand(Problem.N,Problem.D).*beta-Population.decs);
                Ddelta = abs(2*rand(Problem.N,Problem.D).*delta-Population.decs);
                a      = 2*(1-Problem.FE/Problem.maxFE);
                X1     = alpha - (2*a*rand(Problem.N,Problem.D)-a).*Dalpha;
                X2     = beta - (2*a*rand(Problem.N,Problem.D)-a).*Dbeta;
                X3     = delta - (2*a*rand(Problem.N,Problem.D)-a).*Ddelta;
                Population = Problem.Evaluation((X1+X2+X3)/3);
            end
        end
    end
end