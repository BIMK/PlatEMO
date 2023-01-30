classdef RMSProp < ALGORITHM
% <single> <real> <large/none>
% Root mean square propagation
% alpha ---   1 --- Learning rate
% rho   --- 0.9 --- A parameter within [0 1)

%------------------------------- Reference --------------------------------
% T. Tieleman and G. Hinton, Lecture 6.5-rmsprop: Divide the gradient by a
% running average of its recent magnitude, COURSERA: Neural networks for
% machine learning, 2012.
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
            [alpha,rho] = Algorithm.ParameterSet(1,0.9);
            
            %% Generate random solution
            X  = Problem.Initialization(1);
            v0 = zeros(1,Problem.D);

            %% Optimization
            k = 1;
            while Algorithm.NotTerminated(X)
                gk = Problem.CalObjGrad(X.dec);
                v  = rho*v0 + (1-rho)*gk.^2;
                X  = Problem.Evaluation(X.dec-alpha./sqrt(1e-6+v).*gk);
                k  = k + 1;
            end
        end
    end
end