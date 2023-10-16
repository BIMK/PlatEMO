classdef SMOA < ALGORITHM
% <multi> <real> <expensive>
% Supervised multi-objective optimization algorithm
% H --- 2.6e4 --- Size of uniformly distributed L1 unit vector set

%------------------------------- Reference --------------------------------
% T. Takagi, K. Takadama, and H. Sato, Supervised multi-objective
% optimization algorithm using estimation, Proceedings of the IEEE Congress
% on Evolutionary Computation, 2022.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Tomoaki Takagi

% This function requires supervised data like 'DTLZ2_M3_D12.mat' or
% 'DTLZ2_M3_D12.dat'. This function return all evaluated solutions.
% This function do not use random number and evolutionary algorithm.

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            H = Algorithm.ParameterSet(2.6e4);
            
            %% Generate the uniformly distributed L1 unit vector set
            W = UniformPoint(H,Problem.M,'ILD');
            W(W==1e-6) = 0;
            
            %% Generate estimator
            if ~isempty(ver('nnet'))
                estimator = @(W,X,Y) sim(newrbe(X',Y'),W')';
            else
                I = ones(1,M);
                estimator = @(W,X,Y) predictor(W, ...
                    dacefit(X,Y,'regpoly0','corrgauss',I,1e-3*I,1e3*I));
            end
            
            %% Load supervised data
            Population = loadData(Problem);
            Obj = Population.objs;
            Dec = Population.decs;
            
            %% Generate model values
            Y = vecnorm(Obj,1,2); % L1 norm set
            X = Obj./Y;           % L1 unit vector set
            
            %% Search representative vector set
            if Problem.maxFE < length(W)
                Objhat = W.*estimator(W,X,Y);
                Select = subsetSelection(Obj,Objhat,Problem.maxFE);
                W = W(Select,:);
            end
            
            %% Generate decision variables
            decs = zeros(length(W),Problem.D);
            for i = 1 : Problem.D
                decs(:,i) = estimator(W,X,Dec(:,i));
            end
            
            %% Generate population
            Population = [Population,Problem.Evaluation(decs)];
            Algorithm.NotTerminated(Population);
        end
    end
end