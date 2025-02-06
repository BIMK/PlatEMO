function Population = globalReplacement(Population, Offsprings, W, B, T,  PopObj, OffObj)
% Global replacement

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Lucas Farias (email: lucas.farias@unicap.br)

	for i = 1 : length(Offsprings)
		% Calculate the Tchebycheff function value for each weight vector (W) using the offspring's objective values (OffObj)
		tchebycheff_values = max(OffObj(i,:).*W,[],2);		
		
		% Identify the minimum Tchebycheff value, which corresponds to the best solution according to the decomposition method
		best_tchebycheff_value = min(tchebycheff_values);		
		
		% Find the index of the weight vector that corresponds to the best Tchebycheff value
		best_weight_index  = find(tchebycheff_values (:,1) == best_tchebycheff_value);	
		
		% Randomly select a subset of neighborhood solutions from the neighbors of the chosen solution
		P = B(best_weight_index (1),randperm(size(B,2)));					
		
		% Calculate the constraint violation of offspring and P neighborhood
		CVO = sum(max(0,Offsprings(i).con));
		CVP = sum(max(0,Population(P).cons),2);
		
		% Update the solutions in P by Tchebycheff approach
		g_old = max(PopObj(P,:).*W(P,:),[],2);
		g_new = max(repmat(OffObj(i,:),length(P),1).*W(P,:),[],2);
		
		% Replace solutions
		Population(P(find(g_old>=g_new & CVP>=CVO,T))) = Offsprings(i);
	end
end