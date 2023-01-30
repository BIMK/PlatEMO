classdef FRCG < ALGORITHM
% <single> <real> <large/none>
% Fletcher-Reeves conjugate gradient
% beta  --- 0.6 --- A parameter within [0,1] for line search
% sigma --- 0.4 --- A parameter within [0 0.5] for line search

%------------------------------- Reference --------------------------------
% R. Fletcher and C. M. Reeves, Function minimization by conjugate
% gradients, The Computer Journal, 1964, 7(2): 149-154.
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
            k = 0;
            while Algorithm.NotTerminated(X)
                gk = Problem.CalObjGrad(X.dec);
                if mod(k,Problem.D) == 0
                    dk = -gk;
                else
                    betak = (gk*gk')/(g0*g0');
                    dk    = -gk + betak*d0;
                    if gk*dk' >= 0
                        dk = -gk;
                    end
                end
                for m = 0 : 20
                    X1 = Problem.Evaluation(X.dec+beta^m*dk);
                    if X1.obj <= X.obj + sigma*beta^m*gk*dk'
                        break;
                    end
                end
                X  = X1;
                g0 = gk;
                d0 = dk;
                k  = k + 1;
            end
        end
    end
end