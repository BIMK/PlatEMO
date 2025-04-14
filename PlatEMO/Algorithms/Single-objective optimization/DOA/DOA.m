classdef DOA < ALGORITHM
% <2022> <single> <real/integer> <large/none> <constrained/none>
% Dandelion optimization algorithm

%------------------------------- Reference --------------------------------
% S. Zhao, T. Zhang, S. Ma, and M. Chen. Dandelion optimizer: A
% nature-inspired metaheuristic algorithm for engineering applications.
% Engineering Applications of Artificial Intelligence, 2022, 114: 105075.
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
                alpha = rand*((Problem.FE/Problem.maxFE)^2-2*Problem.FE/Problem.maxFE+1);
                % Rise stage
                if randn < 1.5
                    theta    = rand(Problem.N,1)*2*pi-pi;
                    y        = randn(Problem.N,1);
                    InY      = exp(-1/2*log(y).^2)./y./sqrt(2*pi);
                    InY(y<0) = 0;
                    PopDec   = Population.decs + alpha*repmat(cos(theta)./exp(theta),1,Problem.D).*repmat(sin(theta)./exp(theta),1,Problem.D).*repmat(InY,1,Problem.D).*(unifrnd(repmat(Problem.lower,Problem.N,1),repmat(Problem.upper,Problem.N,1))-Population.decs);
                else
                    PopDec = Population.decs.*(1-rand.*((Problem.FE^2-2*Problem.FE+1)/(Problem.maxFE^2-2*Problem.maxFE+1)+1));
                end
                % Decline stage
                beta   = randn(Problem.N,Problem.D);
                PopDec = PopDec - alpha*beta.*(repmat(mean(PopDec,1),Problem.N,1)-alpha*beta.*PopDec);
                % Land stage
                [~,b]  = min(FitnessSingle(Population));
                Elite  = repmat(Population(b).dec,Problem.N,1);
                PopDec = Elite + randn(Problem.N,Problem.D).*(gamma(2.5)*sin(pi*0.75)/gamma(1.25)/1.5/2^0.25)^(1/1.5)./abs(randn(Problem.N,Problem.D)).^(1/1.5)*alpha.*(Elite-PopDec*2*Problem.FE/Problem.maxFE);
                Population = Problem.Evaluation(PopDec);
            end
        end
    end
end