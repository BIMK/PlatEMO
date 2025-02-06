classdef IMCMOEAD < ALGORITHM
% <2024> <multi> <real/integer> <large/none> <constrained/none>
% Inverse modeling constrained MOEA/D
% K --- 10 --- Number of clusters

%------------------------------- Reference --------------------------------
% L. R. C. Farias and A. F. R. Araujo. An inverse modeling constrained
% multi-objective evolutionary algorithm based on decomposition.
% Proceedings of the IEEE International Conference on Systems, Mans and
% Cybernetics, 2024.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Lucas Farias (email: lucas.farias@unicap.br)

    methods
		function main(Algorithm,Problem)
			%% Parameter setting
			K = Algorithm.ParameterSet(10); % Number of clusters
            T = ceil(Problem.N/10); % Size of neighborhood
            
			%% Generate weight vectors
			[W,Problem.N] = UniformPoint(Problem.N,Problem.M);
			
			%% Detect the neighbours of each solution
			B = pdist2(W,W);
			[~,B] = sort(B,2);
			B = B(:,1:T);
			
			%% Generate random population
			Population = Problem.Initialization();			
			idealPoint = min(Population.objs,[],1); % Initial ideal point

			%% Optimization loop
			while Algorithm.NotTerminated(Population)
				% Clustering (k-means)
				[partition,~] = kmeans(Population.objs,K); 
				
				% Modeling and reproduction per cluster
                Offsprings=[];
				for k = unique(partition)'
					parents = Population(partition==k);
                    MatingPool = TournamentSelection(2,size(parents,2),sum(max(0,parents.cons),2));
					Offspring  = Operator(Problem,parents(MatingPool));
					Offsprings = [Offsprings,Offspring];
				end
				
				% Apply constraint handling to Population and Offsprings
                PopObj = applyConstraintHandling(Population, Problem.M);
                OffObj = applyConstraintHandling(Offsprings, Problem.M);
				
				% Update the ideal and nadir points
		        idealPoint = min([idealPoint;OffObj],[],1);
                nadirPoint = max([PopObj;OffObj],[],1);
				
				% Normalize objectives
                PopObj   = (PopObj-repmat(idealPoint,size(PopObj,1),1))./repmat(nadirPoint-idealPoint,size(PopObj,1),1);
                OffObj   = (OffObj-repmat(idealPoint,size(OffObj,1),1))./repmat(nadirPoint-idealPoint,size(OffObj,1),1);
						
				% Performs global replacement based on the Tchebycheff approach
				Population = globalReplacement(Population, Offsprings, W, B, T, PopObj, OffObj);
			end
		end
	end
end