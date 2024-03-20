classdef ISMOA < ALGORITHM
% <multi> <real> <expensive>
% Iterative supervised multi-objective optimization algorithm
% H --- 2.6e4 --- Size of uniformly distributed L1 unit vector set

%------------------------------- Reference --------------------------------
% T. Takagi, K. Takadama, and H. Sato, Pareto front upconvert by iterative
% estimation modeling and solution sampling, Proceedings of The 12th
% Edition of International Conference Series on Evolutionary Multi-
% Criterion Optimization, Lecture Notes in Computer Science, 2023, 13970:
% 218–230.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Tomoaki Takagi

% This function requires solution data like 'DTLZ2_M3_D12_ILD.mat' or
% 'DTLZ2_M3_D12_ILD.dat'. This function return all evaluated solutions.
% This function is deterministic method without random numbers.

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            H = Algorithm.ParameterSet(2.6e4);
            
            %% Generate the uniformly distributed L1 unit vector set
            [W,N] = UniformPoint(H,Problem.M,'ILD');
            W(W==1e-6) = 0;
            
            %% Generate estimator
            if ~isempty(ver('nnet'))
                estimator = @(W,X,Y) sim(newrbe(X',Y'),W')';
            else
                I = ones(1,M);
                estimator = @(W,X,Y) predictor(W, ...
                    dacefit(X,Y,'regpoly0','corrgauss',I,1e-3*I,1e3*I));
            end
            
            %% Load solution data
            Population = loadData(Problem);
            
            %%  Iterative loop
            while size(Population) < Problem.maxFE
                Obj = Population.objs;
                Dec = Population.decs;
                
                % Generate model values
                Y = vecnorm(Obj,1,2); % L1 norm set
                X = Obj./Y;           % L1 unit vector set
                
                % Single solution sampling
                objs = W.*estimator(W,X,Y);
                [~,index] = max(min(pdist2(objs,Obj),[],2));
                
                % Evaluation and upconvert
                dec = zeros(1,Problem.D);
                for i = 1 : Problem.D
                    dec(i) = estimator(W(index,:),X,Dec(:,i));
                end
                dec = min(max(dec,Problem.lower),Problem.upper);
                Population = [Population,Problem.Evaluation(dec)];
            end
            Algorithm.NotTerminated(Population);
        end
    end
end