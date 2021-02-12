classdef BFGS < ALGORITHM
% <single> <real> <large/none>
% A quasi-Newton method proposed by Broyden, Fletcher, Goldfarb, and Shanno
% beta  --- 0.6 --- A parameter within [0,1]
% sigma --- 0.4 --- A parameter within [0 0.5]

%------------------------------- Reference --------------------------------
% D. F. Shanno, Conditioning of quasi-Newton methods for function
% minimization, Mathematics of Computation, 1970, 24: 647-656.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
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
            Bk = eye(Problem.D);
            gk = FiniteDifference(X);
            
            %% Optimization
            while Algorithm.NotTerminated(X)
                dk = -Bk\gk;
                for m = 0 : 20
                    X1 = SOLUTION(X.dec+beta^m*dk');
                    if X1.obj <= X.obj + sigma*beta^m*gk'*dk
                        break;
                    end
                end
                gk1 = FiniteDifference(X1);
                sk  = (X1.dec-X.dec)';
                yk  = gk1 - gk;
                if yk'*sk > 0
                    Bk = Bk - (Bk*sk*sk'*Bk)/(sk'*Bk*sk) + (yk*yk')/(yk'*sk);
                end
                X  = X1;
                gk = gk1;
            end
        end
    end
end