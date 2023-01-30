classdef Adam < ALGORITHM
% <single> <real> <large/none>
% Adaptive moment estimation
% alpha ---     1 --- Learning rate
% beta1 ---   0.9 --- A parameter within [0 1)
% beta2 --- 0.999 --- A parameter within [0 1)

%------------------------------- Reference --------------------------------
% D. P. Kingma and J. Ba, Adam: A method for stochastic optimization, arXiv
% preprint arXiv:1412.6980, 2014.
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
            [alpha,beta1,beta2] = Algorithm.ParameterSet(1,0.9,0.999);
            
            %% Generate random solution
            X  = Problem.Initialization(1);
            m0 = zeros(1,Problem.D);
            v0 = zeros(1,Problem.D);

            %% Optimization
            k = 1;
            while Algorithm.NotTerminated(X)
                gk = Problem.CalObjGrad(X.dec);
                m  = beta1*m0 + (1-beta1)*gk;
                v  = beta2*v0 + (1-beta2)*gk.^2;
                X  = Problem.Evaluation(X.dec-alpha*(m/(1-beta1.^k))./(sqrt(v/(1-beta2.^k))+1e-8));
                k  = k + 1;
            end
        end
    end
end