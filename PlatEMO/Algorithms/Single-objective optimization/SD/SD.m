classdef SD < ALGORITHM
% <single> <real> <large/none>
% Steepest descent
% beta  --- 0.6 --- A parameter within [0,1] for line search
% sigma --- 0.4 --- A parameter within [0 0.5] for line search

%------------------------------- Reference --------------------------------
% S. S. Petrova and A. D. Solov'ev, The origin of the method of steepest
% descent, Historia Mathematica, 1997, 24(4): 361-375.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [beta,sigma] = Algorithm.ParameterSet(0.6,0.4);
            
            %% Generate random solution
            X  = Problem.Initialization(1);

            %% Optimization
            while Algorithm.NotTerminated(X)
                gk = Problem.CalObjGrad(X.dec);
                dk = -gk;
                for m = 0 : 20
                    X1 = Problem.Evaluation(X.dec+beta^m*dk);
                    if X1.obj <= X.obj + sigma*beta^m*gk*dk'
                        break;
                    end
                end
                X = X1;
            end
        end
    end
end