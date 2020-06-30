function PopDec = EGOSelect(Global,Population,L1,L2,Ke,delta,nr)
% Selecting Points for Function Evaluation with the assistance of the
% Gaussian models

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He

    %% Fuzzy clustering the solutions
    [model,centers] = FCMmodel(Global,Population,L1,L2);
    
    %% MOEA/D-DE optimization, where the popsize is set to N, the maximum evaluations is maxeva
	Z       = min(Population.objs,[],1);
    maxeva  = 20*(11*Global.D-1);
    
	%% Generate the weight vectors
	N       = [200 300 595 600 800 800 800 800 800 800];
    [W, N]  = UniformPoint(N,Global.M);
    T       = ceil(N/10);
    B       = pdist2(W,W);
    [~,B]   = sort(B,2);
    B       = B(:,1:T);
    duplicated = randi(length(Population),N,1);
    NewPop  = Population(duplicated);  % Sample N individuals from the initial population
    PopDec  = NewPop.decs; PopObj = NewPop.objs;
    gmin    = inf;
    while (maxeva>0)
       for i = 1 : N
           if rand < delta
               P = B(i,randperm(size(B,2)));
           else
               P = randperm(N);
           end
           OffDec = DE(PopDec(i,:),PopDec(P(1),:),PopDec(P(2),:));
           OffObj = evaluate(Global,OffDec,model,centers);
           Z = min(Z,OffObj);
           g_old = max(abs(PopObj(P,:) - repmat(Z,length(P),1)).*W(P,:),[],2);
           g_new = max(repmat(abs(OffObj-Z),length(P),1).*W(P,:),[],2);
           gmin = min([gmin,min(g_old),min(g_new)]);
           offindex = P(find(g_old>g_new,nr));
           if ~isempty(offindex)
               PopDec(offindex,:) = repmat(OffDec,length(offindex),1); 
               PopObj(offindex,:) = repmat(OffObj,length(offindex),1);
           end
       end      
        maxeva = maxeva - N;
    end

    %% Select the unsimilar candidate solutions
    Q = []; temp = Population.decs;
    for i = 1 : N
        if min(pdist2(real(PopDec(i,:)),real(temp))) > 1e-5
            Q = [Q,i];
            temp = [temp;PopDec(i,:)];
        end
    end
    PopDec = PopDec(Q,:);
    
    %% Kmeans cluster the solutions into Ke clusters and select the solutions with the maximum EI in each cluster
    cindex  = kmeans(real(PopDec),Ke);
    Q = [];
    for i = 1 : Ke
        index = find(cindex == i); 
        temp = PopDec(index,:);
        [tempObj,tempMSE] = evaluate(Global,temp,model,centers);
        K = length(index);
        EI = zeros(K,1);
        for j = 1 : K
            [~,r] = min(max(abs(repmat(tempObj(j,:)-Z,size(W,1),1)).*W,[],2));
            EI(j) = EICal(real(tempObj(j,:)),Z,W(r,:),abs(tempMSE(j,:)),gmin);
        end
        [~,best] = max(EI);
        Q = [Q,index(best)];
    end
    PopDec = PopDec(Q,:);
end

function EI = EICal(Obj,Z,lamda,MSE,Gbest)
% Calculate the expected improvement

    M = size(Obj,2);
    u = lamda.*(Obj - Z);
    sigma2 = lamda.^2.*MSE;
    lamda0 = lamda(1:2); mu0 = u(1:2); sig20 = sigma2(1:2);
	[y,x] = GPcal(lamda0,mu0,abs(sig20));
    if M >= 3
        for i = 3 : M
            lamda0 = [1, lamda(i)]; mu0 = [y,u(i)]; sig20 = [x,sigma2(i)];
            [y,x] = GPcal(lamda0,mu0,abs(sig20));
        end
    end
    EI = (Gbest-y)*normcdf((Gbest-y/sqrt(x))) + sqrt(x)*normpdf((Gbest-y)/sqrt(x));
end

function [y,x] = GPcal(lamda,mu,sig2)
% Calculate the mu (x) and sigma^2 (y) of the aggregation function

    tao = sqrt(lamda(1)^2*sig2(1) + lamda(2)^2*sig2(2));
    alpha = (mu(1)-mu(2))/tao;
    y = mu(1)*normcdf(alpha) + mu(2)*normcdf(-alpha) + tao*normpdf(alpha);
    x = (lamda(1)^2 + sig2(1))*normcdf(alpha) + ...
        (lamda(2)^2 + sig2(2))*normcdf(-alpha) + sum(lamda)*normpdf(alpha);
end

function [PopObj,MSE] = evaluate(Global,PopDec,model,centers)
% Predict the objective vector of the candidate solutions accodring to the
% Euclidean distance from each candidate solution to the evaluated
% solutions
    
    D = pdist2(real(PopDec),real(centers));
    [~,index] = min(D,[],2);
    N = size(PopDec,1);
    PopObj = zeros(N,Global.M);
    MSE = zeros(N,Global.M);
    for i = 1 : N
        for j = 1 : Global.M
            [PopObj(i,j),~,MSE(i,j)] = predictor(PopDec(i,:),model{index(i),j});
        end
    end
end