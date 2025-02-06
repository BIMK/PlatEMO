function new_x = Opt_ETI_FCM(M,D,xlower,xupper,Batch_size,train_x,train_y)
% Maximizing N Subproblems and Selecting Batch of Points 
% Expected Tchebycheff Improvement (ETI)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function was written by Liang Zhao (liazhao5-c@my.cityu.edu.hk).

   %% Fuzzy clustering the solutions
   L1 = 80;
   L2 = 20;
   [GPModels,centers] = GPmodelFCM(train_x,train_y,L1,L2);

   %% Generate the initial weight vectors
    % # of weight vectors：M = 2,  3,  4,  5,  6  
    num_weights = [200,210,295,456,462]; 
    if M <= 3
        [ref_vecs, ~]  = UniformPoint(num_weights(M-1),M); % simplex-lattice design 
    elseif M <= 6
        [ref_vecs, ~]  = UniformPoint(num_weights(M-1),M,'ILD'); % incremental lattice design
    else
        [ref_vecs, ~]  = UniformPoint(500,M); % Two-layered SLD
    end
  
    %% Estimate the Utopian point z
    z       = get_estimation_z(D,xlower,xupper,GPModels,centers,ref_vecs,min(train_y,[],1)); 
    gmin    = get_gmin(train_y,ref_vecs,z); 
   
    %% Using MOEA/D-DE to Maximize ETI
   [pop_ETI,candidate_x,~,~] = MOEAD_ETI(D,xlower,xupper,GPModels,centers,ref_vecs,gmin,z);

    %% Select the unsimilar candidate solutions and build candidate pool Q
    Q = []; Q_ETI = [];  temp = train_x;
    for i = 1 : size(candidate_x,1)
        if min(pdist2(real(candidate_x(i,:)),real(temp))) > 1e-5
            if pop_ETI(i) > 0
                Q = [Q;candidate_x(i,:)]; Q_ETI = [Q_ETI;pop_ETI(i)];
                temp = [temp;candidate_x(i,:)];
            end
        end
    end
 
    new_x = K_means_Batch_Select(Q,Batch_size,candidate_x,Q_ETI) ;
end
 
% >>>>>>>>>>>>>>>>   Auxiliary functions  ==================== 
function gmin = get_gmin(D_objs,ref_vecs,z)
% calculate the minimum of  Tch for each ref_vec
% g(x|w,z) = max{w1(f1-z1),w2(f2-z2),...}
    Objs_translated = D_objs-z; % n*M
    G = ref_vecs(:,1)*Objs_translated(:,1)';  % N*n, f1
    for j = 2:size(ref_vecs,2)
        G = max(G,ref_vecs(:,j)*Objs_translated(:,j)'); % N*n, max(fi,fj)
    end
    gmin = min(G,[],2);  % N*1  one for each weight vector 
end 

function  [pop_ETI,pop_x,pop_mean,pop_std] = MOEAD_ETI(D,xlower,xupper,GPModels,centers,ref_vecs,gmin,z) 
%% using MOEA/D-DE to solve subproblems
    % In order to find the maximum value of ETI for each sub-problem, 
    % it is recommended to set the maximum number of iterations to at least 50.
    %% Parameter setting for MOEA/D-DE
    delta=0.9; nr = 2; maxIter = 50;
    pop_size = size(ref_vecs,1);

    %% neighbourhood   
    T       = ceil(pop_size/10); % size of neighbourhood
    B       = pdist2(ref_vecs,ref_vecs);
    [~,B]   = sort(B,2);
    B       = B(:,1:T);

     %% the initial population for MOEA/D
    pop_x = (xupper-xlower).*lhsdesign(pop_size, D) + xlower; 
    [pop_mean,pop_std] = GPEvaluate_FCM(pop_x,GPModels,centers);
    pop_ETI = get_ETI(pop_mean,pop_std,ref_vecs,gmin,z);   

    for gen = 1 : maxIter-1
       for i = 1 : pop_size
           if rand < delta
               P = B(i,randperm(size(B,2)));
           else
               P = randperm(pop_size);
           end
           %% generate an offspring and calculate its predictive mean and s
           off_x = operator_DE(pop_x(i,:),pop_x(P(1),:),pop_x(P(2),:), xlower,xupper);          
           [off_mean,off_std]= GPEvaluate_FCM(off_x,GPModels,centers);  
            
            ETI_new = get_ETI(repmat(off_mean,length(P),1),repmat(off_std,length(P),1),ref_vecs(P,:),gmin(P),z);
            
            offindex =  find(pop_ETI(P)<ETI_new,nr) ;
            if ~isempty(offindex)
               pop_x(P(offindex),:) = repmat(off_x,length(offindex),1); 
               pop_mean(P(offindex),:) = repmat(off_mean,length(offindex),1);
               pop_std(P(offindex),:) = repmat(off_std,length(offindex),1);
               pop_ETI(P(offindex)) = ETI_new(offindex);
            end
       end      
    end
end

function  ETI = get_ETI(u,sigma,ref_vecs,Gbest,z)
%     g(x|w,z) = max{w1(f1-z1),w2(f2-z2)}  
% calculate the ETI(x|w) at multiple requests, e.g., N  
% u       : N*M  predictive mean
% sigma   : N*M  square root of the predictive variance
% ref_vecs: N*M  weight vectors 
% Gbest   : N*1  
% z       : 1*M  reference point   
    g_mu = ref_vecs.*(u - repmat(z,size(u,1),1));% N*M
    g_sig = ref_vecs.*sigma; % N*M
    % Moment Matching Approximation (MMA)
    g_sig(g_sig<0) = 0; g_sig2 = g_sig.^2; % N*M
     
     % Eq. 18 & Eq. 19 in MOEA/D-EGO
	[mma_mean,mma_sigma2] = app_max_of_2_Gaussian(g_mu(:,1:2),g_sig2(:,1:2)); % f1 & f2
    for i = 3 : size(g_mu,2)
        mu_temp = [mma_mean,g_mu(:,i)]; sig2_temp = [mma_sigma2,g_sig2(:,i)];
        [mma_mean,mma_sigma2] = app_max_of_2_Gaussian(mu_temp,sig2_temp);
    end
    
    mma_std = (sqrt(mma_sigma2));
    Gbest_minus_u = Gbest-mma_mean;
    tau = Gbest_minus_u./mma_std; % n*1

    % Precompute the normal distributions
    normcdf_tau = normcdf(tau);
    normpdf_tau = normpdf(tau);

    ETI = Gbest_minus_u.*normcdf_tau + mma_std.*normpdf_tau;
end

function [u,s] = GPEvaluate_FCM(X,model,centers)
% Predict the objective vector of the candidate solutions accodring to the
% Euclidean distance from each candidate solution to evaluated solutions
    D = pdist2(X,centers);
    [~,index] = min(D,[],2);
    N = size(X,1); % number of samples
    M = size(model,2);% number of objectives
    u = zeros(N,M); % predictive mean
    MSE = zeros(N,M); % predictive MSE
    for i = 1 : N
        for j = 1 : M
            [u(i,j),~,MSE(i,j)] = Predictor(X(i,:),model{index(i),j}); % DACE Kriging toolbox
        end  
    end
     MSE(MSE<0) = 0;
     s = sqrt(MSE);% square root of the predictive variance
end

function [u] = GPEvaluate_mean_FCM(X,model,centers)
% Predict the objective vector of the candidate solutions accodring to the
% Euclidean distance from each candidate solution to evaluated solutions
    D = pdist2(X,centers);
    [~,index] = min(D,[],2);
    N = size(X,1); % number of samples
    M = size(model,2); % number of objectives
    u = zeros(N,M); % predictive mean
    for i = 1 : N
        for j = 1 : M
            [u(i,j)] = Predictor(X(i,:),model{index(i),j}); % DACE Kriging toolbox
        end
    end
end

function [y,s2] = app_max_of_2_Gaussian(mu,sig2)
% Calculate  Eq. 18 & Eq. 19 in MOEA/D-EGO 
% n requests
% mu is N*2
% sig2 is N*2
    tao = sqrt(sum(sig2,2));  % N*1
    alpha = (mu(:,1)-mu(:,2))./tao;  % N*1
    % Eq. 16 / Eq. 18
    y = mu(:,1).*normcdf(alpha) + mu(:,2).*normcdf(-alpha) + tao.*normpdf(alpha);  % N*1
    % There is a typo in Eq. 17.  See Appendix B of MOEA/D-EGO.
    % It should be $$ +(\mu_1+\mu_2) \tau \varphi(\alpha)$$
    s2 = (mu(:,1).^2 + sig2(:,1)).*normcdf(alpha) + ...
        (mu(:,2).^2 + sig2(:,2)).*normcdf(-alpha) + (sum(mu,2)).*tao.*normpdf(alpha); 
%     s2 = (mu(:,1).^2 + sig2(:,1)).*normcdf(alpha) + ...
%         (mu(:,2).^2 + sig2(:,2)).*normcdf(-alpha) + (sum(mu,2)).*normpdf(alpha); 
    s2 = s2 - y.^2;
    s2(s2<0) = 0;
end

function  z = get_estimation_z(D, xlower,xupper,GPModels,centers,ref_vecs,z) 
% min(\mu_1(x),...,\mu_m(x))^T
% Utilize MOEA/D to minimize the GP posterior mean and determine the utopian
% point during optimization. Alternatively, other multi-objective optimization 
% algorithms such as NSGA-II can also be employed.

    delta=0.9; nr = 2; 
    maxIter = 100; 
    pop_size = size(ref_vecs,1);

    %% neighbourhood   
    T       = ceil(pop_size/10); % size of neighbourhood
    B       = pdist2(ref_vecs,ref_vecs);
    [~,B]   = sort(B,2);
    B       = B(:,1:T);

    %% the initial population for MOEA/D
    pop_x = (xupper-xlower).*lhsdesign(pop_size, D) + xlower;
    pop_mean = GPEvaluate_mean_FCM(pop_x,GPModels,centers);
    z       = min(min(pop_mean,[],1),z);
    for gen = 1 : maxIter-1
       for i = 1 : pop_size
           if rand < delta
               P = B(i,randperm(size(B,2)));
           else
               P = randperm(pop_size);
           end
           %% generate an offspring and calculate its predictive mean and s
           off_x = operator_DE(pop_x(i,:),pop_x(P(1),:),pop_x(P(2),:), xlower,xupper);    
           [off_mean]= GPEvaluate_mean_FCM(off_x,GPModels,centers);  
           
            z = min(z,off_mean);        
            g_old = max((pop_mean(P,:) - repmat(z,length(P),1)).*ref_vecs(P,:),[],2);
            g_new = max(repmat((off_mean-z),length(P),1).*ref_vecs(P,:),[],2);
             
           % Update the solutions in P
           offindex = P(find(g_old>g_new,nr));
           if ~isempty(offindex)
               pop_x(offindex,:) = repmat(off_x,length(offindex),1); 
               pop_mean(offindex,:) = repmat(off_mean,length(offindex),1);
           end
       end      
    end
end
function SelectDecs = K_means_Batch_Select(Q,Batch_size,candidate_x,Q_ETI) 
     batch_size = min(Batch_size,size(Q,1));% in case Q is smaller than Batch size
    
    if batch_size == 0
        Qb = randperm(size(candidate_x,1),Batch_size);
        SelectDecs = candidate_x(Qb,:);
    else
        cindex  = kmeans(Q,batch_size);
        Qb = [];
        for i = 1 : batch_size
            index = find(cindex == i); 
            [~,best] = max(Q_ETI(index));
            Qb = [Qb,index(best)];
        end
        SelectDecs = Q(Qb,:);
    end
    % when Q is smaller than batch size
    if size(SelectDecs,1) < Batch_size
        Qb = randperm(size(candidate_x,1), Batch_size - size(SelectDecs,1));
        SelectDecs = [SelectDecs;candidate_x(Qb,:)];
    end
end
% >>>>>>>>>>>>>>>>    functions in PlatEMO ====================
function Offspring = operator_DE(Parent1,Parent2,Parent3, xlower,xupper)
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