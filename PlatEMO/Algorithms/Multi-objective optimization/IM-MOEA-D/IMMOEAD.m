classdef IMMOEAD < ALGORITHM
% <multi> <real/integer> <large/none>
% Inverse modeling multiobjective evolutionary algorithm based on
% decomposition
% K --- 10 --- Number of reference vectors

%------------------------------- Reference --------------------------------
% L. R. C. Farias and A. F. R. Araujo, IM-MOEA/D: An inverse modeling
% multi-objective evolutionary algorithm based on decomposition,
% Proceedings of the IEEE International Conference on Systems, Mans and
% Cybernetics, 2021.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Lucas Farias

    methods
		function main(Algorithm,Problem)
			%% Parameter setting
			K = Algorithm.ParameterSet(10);
            T = ceil(Problem.N/10); % Size of neighborhood
            
			%% Generate weight vectors
			[W,Problem.N] = UniformPoint(Problem.N,Problem.M);
			
			%% Detect the neighbours of each solution
			B = pdist2(W,W);
			[~,B] = sort(B,2);
			B = B(:,1:T);
			
			%% Generate random population
			Population    = Problem.Initialization();
			Z = min(Population.objs,[],1);
			
			%% Optimization
			while Algorithm.NotTerminated(Population)
				[partition,~] = kmeans(Population.objs,K); 
				Offsprings=[];
				% Modeling and reproduction
				for k = unique(partition)'
					Offspring  = Operator(Problem,Population(partition==k));
					Offsprings = [Offsprings,Offspring];
				end
				
				% Update the ideal point
				Z = min(Z,min(Offsprings.objs));		
				
				for i = 1 : length(Offsprings)
					% Global Replacement
					all_g_TCH=max(abs((Offsprings(i).obj-repmat(Z,Problem.N,1)).*W),[],2);
					best_g_TCH=min(all_g_TCH);
					Chosen_one = find(all_g_TCH(:,1)==best_g_TCH);
					P = B(Chosen_one(1),randperm(size(B,2)));
					
					% Update the solutions in P by Tchebycheff approach
					g_old = max(abs(Population(P).objs-repmat(Z,length(P),1)).*W(P,:),[],2);
					g_new = max(repmat(abs(Offsprings(i).obj-Z),length(P),1).*W(P,:),[],2);
					Population(P(find(g_old>=g_new,T))) = Offsprings(i);	
				end
			end
		end
	end
end