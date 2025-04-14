classdef WOA < ALGORITHM
% <2016> <single> <real/integer> <large/none> <constrained/none>
% Whale optimization algorithm

%------------------------------- Reference --------------------------------
% S. Mirjalili and A. Lewis. The whale optimization algorithm. Advances in
% Engineering Software, 2016, 95: 51-67.
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
                A = 2*(2-2*Problem.FE/Problem.maxFE)*rand(Problem.N,1) - (2-2*Problem.FE/Problem.maxFE);
                C = 2*rand(Problem.N,1);
                type1  = rand(Problem.N,1) >= 0.5;
                type2  = ~type1 & abs(A)<1;
                type3  = ~type1 & abs(A)>=1;
                [~,b]  = min(FitnessSingle(Population));
                PopDec = Population.decs;
                % Spiral updating position
                l = repmat(rand(sum(type1),1)*2-1,1,Problem.D);
                PopDec(type1,:) = abs(repmat(PopDec(b,:),sum(type1),1)-PopDec(type1,:)).*exp(l).*cos(2*pi*l) + repmat(PopDec(b,:),sum(type1),1);
                % Encircling prey
                PopDec(type2,:) = repmat(PopDec(b,:),sum(type2),1) - repmat(A(type2),1,Problem.D).*abs(C(type2)*PopDec(b,:)-PopDec(type2,:));
                % Search for prey
                r = randperm(Problem.N,sum(type3));
                PopDec(type3,:) = PopDec(r,:) - repmat(A(type3),1,Problem.D).*abs(repmat(C(type3),1,Problem.D).*PopDec(r,:)-PopDec(type3,:));
                Population = Problem.Evaluation(PopDec);
            end
        end
    end
end