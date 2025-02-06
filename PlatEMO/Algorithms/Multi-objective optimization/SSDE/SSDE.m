classdef SSDE < ALGORITHM
% <2024> <multi/many> <real/integer> <constrained/none> <expensive>
% Self-organized surrogate-assisted differential evolution
% num_nodes ---     --- Number of neurons in each dimension of the latent space
% eta0      --- 0.2 --- Initial learning rate
% sigma0    ---     --- Size of neighborhood mating pools

%------------------------------- Reference --------------------------------
% A. F. R. Araújo, L. R. C. Farias, and A. R. C. Gonçalves. Self-organizing
% surrogate-assisted non-dominated sorting differential evolution. Swarm
% and Evolutionary Computation, 2024, 91: 101703.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Lucas Farias (email: lrcf@cin.ufpe.br)

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [num_nodes,eta0,sigma0] = Algorithm.ParameterSet(Problem.N, 0.2, Problem.N);

            %% Generate random population
			Population = Problem.Initialization();

			%% Initialize the sample set
			Samples = Population;
			
			%% Initialize the SOM weight vectors
            a = 0.001;
            b = 0.5;
            W = a.*randn(num_nodes,Problem.D+Problem.M) + b;
			winning_weights = false(1,Problem.N);			
			
			%% Position of each neuron and neighborhood
			D = arrayfun(@(S)1:S,floor(num_nodes),'UniformOutput',false);
			eval(sprintf('[%s]=ndgrid(D{:});',sprintf('c%d,',1:length(D))))
			eval(sprintf('V=[%s];',sprintf('c%d(:),',1:length(D))))			
			LDis = pdist2(V,V); % Distance between each two neurons in latent space
			
            %% Optimization
            while Algorithm.NotTerminated(Population)
				%% Training
                if size(Samples,2) >= Problem.N
                    W = Training(Problem, W, LDis, Samples, num_nodes, eta0, sigma0, winning_weights);
                    winning_weights = false(1,num_nodes);
					Samples = [];
                end
                
                %% Evolutionary process
                [Population, new_samples, winning_weights] = Operator(Problem, Population, W, winning_weights);
				Samples = [Samples, new_samples];
            end
        end
    end
end