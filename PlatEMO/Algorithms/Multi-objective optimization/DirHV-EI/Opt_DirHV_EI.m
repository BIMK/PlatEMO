function new_x = Opt_DirHV_EI(M,D,xlower,xupper,GPModels,train_y_nds,Batch_size)
% Maximizing N Subproblems and Selecting Batch of Query Points 
% Expected Direction-based Hypervolume Improvement (DirHV-EI, denoted as EI_D)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function was written by Liang Zhao.
% https://github.com/mobo-d/DirHV-EGO

   %% Generate the initial reference vectors
    % # of weight vectors：M = 2,  3,  4,  5,  6  
    num_weights = [200,210,295,456,462]; 
    if M <= 3
        [ref_vecs, ~] = UniformPoint(num_weights(M-1), M); % simplex-lattice design 
    elseif M <= 6
        [ref_vecs, ~] = UniformPoint(num_weights(M-1), M,'ILD'); % incremental lattice design
    else
        [ref_vecs, ~] = UniformPoint(500, M);   % Two-layered SLD
    end
 
   %% Utopian point
    z = -0.01.*ones([1, M]); % Adaptively adjusting Z may lead to better performance.

   %% Calculate the Intersection points and Direction vectors
    [xis,dir_vecs] = get_xis(train_y_nds,ref_vecs,z);  
  
   %% Use MOEA/D-GR to maximize DirHV-EI over a set of direction vectors in a collaborative manner
    [~,candidate_x,candidate_mean,candidata_std] = MOEAD_GR_(D,xlower,xupper,GPModels,dir_vecs,xis);
 
    %% discard the duplicate candidates
    [candidate_x,ia,~] = unique(candidate_x,'rows'); 
    candidate_mean     = candidate_mean(ia,:); candidata_std = candidata_std(ia,:);
 
    %% Compute EI_D for all the points in Q
    DirHVEIs = zeros(size(candidate_x,1),size(dir_vecs,1));
    for j = 1 : size(candidate_x,1)
	    DirHVEIs(j,:) = get_DirHVEI(repmat(candidate_mean(j,:),size(dir_vecs,1),1),repmat(candidata_std(j,:),size(dir_vecs,1),1),xis); 
    end
    %% find q solutions with the greedy algorithm
    Qb    = subset_selection(DirHVEIs,Batch_size);  
    new_x = candidate_x(Qb,:); 
end

% >>>>>>>>>>>>>>>>   Algorithm 2 & Algorithm 3  ====================
function [pop_dirhvei,pop_x,pop_mean,pop_std] = MOEAD_GR_(D,xlower,xupper,GPModels,dir_vecs,xis)
%% Algorithm 2: using MOEA/D-GR to solve subproblems
    % In order to find the maximum value of DirHV-EI for each sub-problem, 
    % it is recommended to set the maximum number of iterations to at least 50.
    maxIter  = 50; 
    pop_size = size(dir_vecs,1);
    %% neighbourhood   
    T     = ceil(pop_size/10); % size of neighbourhood: 0.1*N
    B     = pdist2(dir_vecs,dir_vecs);
    [~,B] = sort(B,2);
    B     = B(:,1:T);
	
    % the initial population for MOEA/D
    pop_x = (xupper-xlower).*lhsdesign(pop_size, D) + xlower; 
    [pop_mean,pop_std] = GPEvaluate(pop_x,GPModels); % calculate the predictive means and variances
    pop_dirhvei = get_DirHVEI(pop_mean,pop_std,xis); 
	
	% optimization
    for gen = 1 : maxIter-1
       for i = 1 : pop_size    
           if rand < 0.8 % delta
               P = B(i,randperm(size(B,2)));
           else
               P = randperm(pop_size);
           end
           %% generate an offspring and calculate its predictive mean and variance 
           off_x = operator_DE(pop_x(i,:),pop_x(P(1),:),pop_x(P(2),:), xlower,xupper); 
           [off_mean,off_std] = GPEvaluate(off_x,GPModels);  
            
           %% Global Replacement  MOEA/D-GR
           % Find the most appropriate subproblem and its neighborhood
            DirHVEIs = get_DirHVEI(repmat(off_mean,pop_size,1),repmat(off_std,pop_size,1),xis);
            [~,best_index] = max(DirHVEIs);

            P = B(best_index,:); % replacement neighborhood
            % Update the solutions in P
            offindex = P(pop_dirhvei(P)<DirHVEIs(P));
           if ~isempty(offindex)
               pop_x(offindex,:)     = repmat(off_x,length(offindex),1); % pop_x: N*D
               pop_mean(offindex,:)  = repmat(off_mean,length(offindex),1); % pop_mean: N*M
               pop_std(offindex,:)   = repmat(off_std,length(offindex),1); % pop_std: N*M
               pop_dirhvei(offindex) = DirHVEIs(offindex); % pop_dirhvei: N*1
           end
       end   
    end
end

function Qb = subset_selection(DirHVEIs,Batch_size)
%% Algorithm 3: Submodularity-based Batch Selection
    [L,N] = size(DirHVEIs);
    Qb    = [];
    temp  = DirHVEIs;
    beta  = zeros([1,N]); 
    for i = 1 : Batch_size
        [~,index] = max(sum(temp,2));
        Qb   = [Qb,index];
        beta = beta + temp(index,:);
        % temp: [EI_D(x|\lambda) - beta]_+
        temp = DirHVEIs-repmat(beta,L,1);
        temp(temp < 0) = 0;   
    end
end

% >>>>>>>>>>>>>>>>   Auxiliary functions  ====================
function [xis,dir_vecs] = get_xis(train_y_nds,ref_vecs,z)  
% % calculate the minimum of mTch for each direction
% % L(A,w,z) = min max (A-z)./Lambda % used for Eq.11 
% % In:
% % train_y_nds  : n*M  observed objectives 
% % ref_vecs     : N*M  weight vectors
% % z            : 1*M  reference point
% % Return:
% % xis          : N*1  intersection points, one for each direction vector 
% % dir_vecs     : N*M  direction vectors 

    [~,M] = size(ref_vecs); % N: # of subproblems, M: # of objectives
    %% Eq. 25, the L2 norm of each direction vector should be 1
    temp      = 1.1.*ref_vecs - z;
    norm_temp = sqrt(sum(temp.^2, 2));
    dir_vecs  = temp ./ norm_temp;
    
	%% Eq. 11, compute the intersection points
    div_dir = 1./dir_vecs; % N*M
    train_y_translated = train_y_nds-z; % n*M
    G = div_dir(:,1)*train_y_translated(:,1)'; % N*n, f1
    for j = 2 : M
        G = max(G,div_dir(:,j)*train_y_translated(:,j)'); % N*n, max(fi,fj)
    end
    % minimum of mTch for each direction vector
    Lmin = min(G,[],2); % N*1  one for each direction vector 
    % N*M  Intersection points
    xis = z + Lmin.*dir_vecs;% Eq.11
end 

function DirHVEI = get_DirHVEI(u,sigma,xis)
% calculate the EI_D(x|lambda) at multiple requests 
% u     : N*M  predictive mean
% sigma : N*M  square root of the predictive variance
% Xis   : N*M  Intersection points 
% % Eq. 23 in Proposition 5
    xi_minus_u = xis-u; % N*M
    tau = xi_minus_u./sigma;  % N*M   

    % Precompute the normal distributions
    normcdf_tau = normcdf(tau);
    normpdf_tau = normpdf(tau);

    temp = xi_minus_u .* normcdf_tau + sigma .* normpdf_tau; % N*M
    DirHVEI = prod(temp,2);
end 
 
function [u,s] = GPEvaluate(X,model)
% Predict the GP posterior mean and std at a set of the candidate solutions 
    N   = size(X,1);        % number of samples
    M   = length(model);    % number of objectives
    u   = zeros(N,M);       % predictive mean
    MSE = zeros(N,M);       % predictive MSE
    if N == 1 
        for j = 1 : M
            [u(:,j),~,MSE(:,j)] = Predictor(X,model{1,j}); % DACE Kriging toolbox
        end
        MSE(MSE<0) = 0;
    else
        for j = 1 : M
            [u(:,j),MSE(:,j)] = Predictor(X,model{1,j}); % DACE Kriging toolbox
        end
        MSE(MSE<0) = 0;
    end
   s = sqrt(MSE);% square root of the predictive variance
end

% >>>>>>>>>>>>>>>>    functions in PlatEMO ====================
function Offspring = operator_DE(Parent1,Parent2,Parent3,xlower,xupper)
%OperatorDE - The operator of differential evolution.

    %% Parameter setting
    [CR,F,proM,disM] = deal(1,0.5,1,20);
    [N,D] = size(Parent1);

    %% Differental evolution
    Site = rand(N,D) < CR;
    Offspring       = Parent1;
    Offspring(Site) = Offspring(Site) + F*(Parent2(Site)-Parent3(Site));

    %% Polynomial mutation
    Lower = repmat(xlower,N,1);
    Upper = repmat(xupper,N,1);
    Site  = rand(N,D) < proM/D;
    mu    = rand(N,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5; 
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
end