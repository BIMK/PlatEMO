function W = Training(Problem, W, LDis, Samples, num_nodes, eta0, sigma0, winning_weights)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Lucas Farias (email: lrcf@cin.ufpe.br)

    %% Memory-based reinitialization
	if any(winning_weights) == true                    
        %% If there is a node that did not win
        if sum(winning_weights) < num_nodes   
            %% Select indexes based on diversity criteria
            [FrontNo,~] = NDSort(W(winning_weights,Problem.D+1:end),sum(winning_weights));

			%% Calculate the crowding distance of each solution
            CrowdDis   = CrowdingDistance(W(winning_weights,Problem.D+1:end),FrontNo); 
            [~,p]      = sort(CrowdDis);
            factor     = 1 : length(CrowdDis);
            factor(p)  = factor; % each winning node receives a weight according to the diversity criterion
            chosen_idx = randsample(sum(winning_weights),2*(num_nodes-sum(winning_weights)),true,factor);
            chosen_idx = reshape(chosen_idx,num_nodes-sum(winning_weights),2);

            % linear combination: calculate average and add noise to new nodes
            W(~winning_weights,:) = ((W(chosen_idx(:,1),:) + W(chosen_idx(:,2),:) ) ./ 2) + 0.001.*randn(num_nodes-sum(winning_weights),Problem.D+Problem.M);

            % repair if it exceeds the upper limit, 1, or lower limit, 0
            W(:,1:Problem.D) = max(min(W(:,1:Problem.D),1),0);
        end
	end
		
	%% Normalize Sample Set
	Samples = [Samples.decs,Samples.objs];
    Samples(:,1:Problem.D) = rescale(Samples(:,1:Problem.D),'InputMin',Problem.lower,'InputMax',Problem.upper);

    % reset parameters
	win_count_set = zeros(1,num_nodes);

    % SOM surrogate model training
	for epoch = 1 : 50
		% shuffle the samples
		randpos = randperm(size(Samples,1));        
		% Neighborhood radius 
		sigma = sigma0*exp((-win_count_set/size(Samples,1))); 
		% Learning rate		
		eta = eta0*exp((-win_count_set/size(Samples,1))); 
		
        for i = 1 : size(Samples,1) 
			s = randpos(i);
            
            % define winning node, u1, from input sample, s, using quadratic Euclidean distance
			[~,u1] = min(pdist2(Samples(s,1:end-Problem.M),W(:,1:end-Problem.M)));

			% If it is the node's first win, assign the sample to the node
            if win_count_set(u1) == 0
				W(u1,:) = Samples(s,:);
            end
            
            % limit node win counter to current epoch
            if win_count_set(u1) < epoch
                win_count_set(u1) = win_count_set(u1) + 1;
            end
            
            % update winning node and its neighborhood
			U      = LDis(u1,:) < sigma;
			W(U,:) = W(U,:) + repmat(eta(U)',1,size(W,2)).*repmat(exp(-LDis(u1,U))',1,size(W,2)).*(repmat(Samples(s,:),sum(U),1)-W(U,:));
        end
	end
end