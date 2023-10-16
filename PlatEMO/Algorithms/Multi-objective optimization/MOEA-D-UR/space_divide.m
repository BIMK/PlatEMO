function [Population,Z]=space_divide(Problem,EP,Population,Z,W_URP,W,K,mini_generation)
% Division of the Objective Space

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Lucas Farias

	[idx,~] = kmeans(Population.objs,K); 
	
	allInd = [Population EP];

	%% For each group
	for current_group_index=1:max(idx)		
		%% Define the current group
		currentGroup=Population(idx(:)==current_group_index);
		group_size=length(currentGroup);
        if length(unique(currentGroup)) > 2
            %% Generate the normalized weight vectors
            W_normalized=W_URP(idx(:)==current_group_index,:);
            w_maximum=max(W_normalized);
            w_minimums=min(W_normalized);
            W_normalized= (W_URP.*repmat((w_maximum-w_minimums),length(W_URP),1)) + repmat(w_minimums,length(W),1); %normalization in variable range(x,y)

            %% Associate individuals to weights
            currentGroup=associate_newW2allInd(allInd,W_normalized,Z);

            [candidates,Z]=mini_gen(Problem,currentGroup,Z,group_size,mini_generation);
            allInd=[allInd candidates];
        end
	end
	Population = combine_allInd2Population(Problem,Population,allInd,W,Z);
end

function Population = associate_newW2allInd(candidates,W,Z)   
	Combine = candidates;
    CombineObj = abs(Combine.objs-repmat(Z,length(Combine),1));
    g = zeros(length(Combine),size(W,1));
    for i = 1 : size(W,1)
        g(:,i) = max(CombineObj.*repmat(W(i,:),length(Combine),1),[],2);
    end
    % Choose the best solution for each subproblem
    [~,best]   = min(g,[],1);
    Population = Combine(best);
end

function [currentGroup,Z] = mini_gen(Problem,currentGroup,Z,group_size,mini_generation)
	generation = 1;
	while generation <= mini_generation
		
		% Choose the parents
		parents=unique(currentGroup);
		MatingPool = randi(length(parents),1,group_size);
		% Generate offsprings
        Offsprings  = OperatorGA(Problem,parents(MatingPool));
        
		% Update the ideal point
		Z = min([Z;Offsprings.objs],[],1);
		
		currentGroup=[currentGroup Offsprings];
		
        generation=generation+1;
	end
end

function Population = combine_allInd2Population(Problem,Population,candidates,W,Z)
   
    Combine=unique(candidates);
    if length(Combine) < Problem.N
        Combine = candidates;
    end
    CombineObj = abs(Combine.objs-repmat(Z,length(Combine),1));
    g = zeros(length(Combine),size(W,1));
    for i = 1 : size(W,1)
        g(:,i) = max(CombineObj.*repmat(W(i,:),length(Combine),1),[],2);
    end
    % Choose the best solution for each subproblem
    [~,best]   = min(g,[],1);
    candidates = Combine(best);
    candidates=unique(candidates);
	
    for i=1:length(candidates)        
        % Global Replacement
        all_g_TCH=max(abs((candidates(i).obj-repmat(Z,Problem.N,1)).*W),[],2);
        best_g_TCH=min(all_g_TCH);
        Chosen_one = find(all_g_TCH(:,1)==best_g_TCH);            
        % Update the solutions in P by Tchebycheff approach        
        if Population(Chosen_one) ~= candidates(i)
            g_old = max(abs(Population(Chosen_one).objs-repmat(Z,length(Chosen_one),1)).*W(Chosen_one,:),[],2);
            g_new = max(repmat(abs(candidates(i).obj-Z),length(Chosen_one),1).*W(Chosen_one,:),[],2);
            Population(Chosen_one(g_old>=g_new)) = candidates(i);
        end
    end
end