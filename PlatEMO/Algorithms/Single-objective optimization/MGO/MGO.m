classdef MGO < ALGORITHM
% <2022> <single> <real/integer> <large/none> <constrained/none>
% Mountain gazelle optimizer

%------------------------------- Reference --------------------------------
% B. Abdollahzadeh, F. S. Gharehchopogh, N. Khodadadi, and S. Mirjalili.
% Mountain gazelle optimizer: A new nature-inspired metaheuristic algorithm
% for global optimization problems. Advances in Engineering Software, 2022,
% 174: 103282.
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
            Randi      = @()repmat(randi(2,Problem.N,1),1,Problem.D);
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                Cof(1,:) = repmat(-Problem.FE/Problem.maxFE+rand,1,Problem.D);
                Cof(2,:) = (-1-Problem.FE/Problem.maxFE).*randn(1,Problem.D);
                Cof(3,:) = rand(1,Problem.D);
                N3       = randn(1,Problem.D);
                Cof(4,:) = N3.*randn(1,Problem.D).^2.*cos(rand*2*N3);
                for i = 1 : Problem.N
                    BH(i,:) = mean(Population(randperm(Problem.N,ceil(Problem.N/3))).decs,1);
                end
                TSM = repmat(Population(1).dec,Problem.N,1) - abs((Randi().*BH-Randi().*Population.decs).*randn(Problem.N,Problem.D)*exp(2-2*Problem.FE/Problem.maxFE)).*Cof(randi(end,1,Problem.N),:);
                MH  = BH + Cof(randi(end,1,Problem.N),:) + (randi(2,Problem.N,1)*Population(1).dec-Randi().*Population(randi(end,1,Problem.N)).decs).*Cof(randi(end,1,Problem.N),:);
                BMH = Population.decs - (abs(Population.decs)+repmat(abs(Population(1).dec),Problem.N,1)).*(2*repmat(rand(Problem.N,1),1,Problem.D)-1) + (randi(2,Problem.N,1)*Population(1).dec-Randi().*BH).*Cof(randi(end,1,Problem.N),:);
                MSF = unifrnd(repmat(Problem.lower,Problem.N,1),repmat(Problem.upper,Problem.N,1));
                Population = [Population,Problem.Evaluation([TSM;MH;BMH;MSF])];
                [~,rank]   = sort(FitnessSingle(Population));
                Population = Population(rank(1:Problem.N));
            end
        end
    end
end