function MESMO(Global)
% <algorithm> <GA> <expensive!>
% Max-value Entropy Search for Multi-Objective Bayesian Optimization
% nK             --- 10 --- Number of samples to approximate the MES
% learn_interval ---  5 --- Iterations to update hyperparameters 

%------------------------------- Reference --------------------------------
% S. Belakaria, A. Deshwal, J. R. Doppa, Max-value Entropy Search for
% Multi-Objective Bayesian Optimization, Proceedings of the 33rd Conference
% on Neural Information Processing Systems, 2019, 7825-7835.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Youwei He (1554748356@qq.com)
	
	[nK,learn_interval] = Global.ParameterSet(10,5);
    %% Generate the initial design points
    % number of design variables
    D = Global.D;
    % number of objective functions
    M = Global.M;
    % number of initial design points
    N = Global.N;%11*Global.D-1
    % generate initial design points using Latin Hypercube sampling
    PopDec = repmat(Global.upper-Global.lower,N,1).*lhsamp(N,D)+repmat(Global.lower,N,1);
    % calculate initial design points
    Population   = INDIVIDUAL(PopDec);
    % parameters;
    nM = 1;  % must not be modified
    nFeatures = 1000;  iter = 0;

    %% Optimization
    while Global.NotTermination(Population)
        PopDec = Population.decs;
        PopObj = Population.objs;%% objectives to be maximized
        % sample from the posterior distribution of the hyper-parameters     
        if ~mod(iter,learn_interval) 
        l = cell(M, 1); sigma = cell(M, 1); sigma0 = cell(M, 1); KernelMatrixInv = cell(M, nM);
            for i=1:M
                [ l{i}, sigma{i}, sigma0{i} ] = sampleHypers(PopDec, PopObj(:,i), nM, []);
            end
        end
        for i=1:M
            for j = 1 : nM
                lKmm = l{i}(j,:)'; sigmaKmm = sigma{i}(j); sigma0Kmm = sigma0{i}(j);
                KernelMatrix = computeKmm(PopDec, lKmm, sigmaKmm, sigma0Kmm);
                KernelMatrixInv{ i,j } = chol2invchol(KernelMatrix);
            end
        end
        % get the infill point
        PopDec = MESMO_Choose(nM, nK, PopDec, PopObj, KernelMatrixInv,  ...
            sigma0, sigma, l, Global.lower, Global.upper, nFeatures);
        Population = [Population,INDIVIDUAL(PopDec)];
        iter = iter+1;
    end
end

%% Whole procedure for choosing the next point
function Best_x = MESMO_Choose(nM, nK, xx, yy, KernelMatrixInv, sigma0, sigma, l, xmin, xmax, nFeatures)
D = size(xx,2);
nObj = size(yy,2);
% sample the objective functions and get the PFs using cheap
% Multi-objective optimization method like NSGA-II(RVEA)
PFs = GetPFs(nM, nK, xx, yy, sigma0, sigma, l, xmin, xmax, nFeatures, nObj);
% Define the acquisition function of MOMES.
acfun = @(x) evaluateMOMES(x, PFs, xx, yy, KernelMatrixInv, l, sigma, sigma0);
% Optimize the acquisition function using GA
% ---------  by GA                     --------
options = optimoptions('ga','MaxGenerations',400,'MaxStallGenerations',50, ...
        'PopulationSize',min(40*D,120),'CrossoverFraction',0.8,'MigrationFraction',0.2, ...
        'Display','off','UseVectorized',true);
Best_x = ga(acfun, D, [],[],[],[], xmin, xmax, [], options);  
end

%% Get Pfs using NSGA-II
function PFs = GetPFs(nM, nK, xx, yy, sigma0Obj, sigmaObj, lObj, xmin, xmax, nFeatures, nObj)
% nM is the number of sampled GP hyper-parameter settings.
% nK is the number of sampled maximum values.
% xx, yy are the current observations.
% sigma0, sigma, l are the hyper-parameters of the Gaussian kernel.
% xmin, xmax are the lower and upper bounds for the search space.
% nFeatures is the number of random features sampled to approximate the GP
D = size(xx,2);
PFs = cell(nM,nK);
for i=1:nM
    for j=1:nK
        l = zeros(nObj,D); sigma0 = zeros(nObj,1); sigma = zeros(nObj,1);
        for k=1: nObj
            sigma0(k) = sigma0Obj{k}(i);
            sigma(k) = sigmaObj{k}(i);
            l(k,:) = lObj{k}(i,:);
        end
        ObjFuns = @(x)sampledObjFuns(x, xx, yy, sigma0, sigma, l, nFeatures, nObj);
        % NSGS-II optimization
        objNSGA2 = NSGA_2( ObjFuns, xmin,xmax, D, 'n_pop',20*D,'max_gen',300);
        PFs{i,j} = objNSGA2.y;
    end
end
end

%% sample objective functions
function [obj,Cons] = sampledObjFuns(x, xx, yy, sigma0, sigma, l, nFeatures, nObj)
d = size(xx,2);
nx = size(x,1);
obj = zeros(nx,nObj);
for k=1:nObj
        % Draw weights for the random features.
        W = randn(nFeatures, d) .* repmat(sqrt(l(k,:)), nFeatures, 1);% nFeatures by d
        b = 2 * pi * rand(nFeatures, 1);% nFeatures by 1

        % Compute the features for xx. % nFeatures by nx
        Z = sqrt(2 * sigma(k) / nFeatures) * cos(W * xx' + repmat(b, 1, size(xx, 1)));

        % Draw the coefficient theta.
        noise = randn(nFeatures, 1);
        if (size(xx, 1) < nFeatures)
            % We adopt the formula $theta \sim \N(Z(Z'Z + \sigma^2 I)^{-1} y, 
            % I-Z(Z'Z + \sigma^2 I)Z')$.
            Sigma = Z' * Z + sigma0(k) * eye(size(xx, 1));% nx by nx
            mu = Z*chol2invchol(Sigma)*yy(:,k);% nFeatures by 1
            [U, D] = eig(Sigma);
            D = diag(D);% nx by 1
            R = (sqrt(D) .* (sqrt(D) + sqrt(sigma0(k)))).^-1;% nx by 1
            theta = noise - (Z * (U * (R .* (U' * (Z' * noise))))) + mu;% nFeatures by 1
        else
            % $theta \sim \N((ZZ'/\sigma^2 + I)^{-1} Z y / \sigma^2,
            % (ZZ'/\sigma^2 + I)^{-1})$.---incorrect
            Sigma = chol2invchol(Z*Z' / sigma0(k) + eye(nFeatures));% nFeatures by nFeatures
            mu = Sigma * Z * yy(:,k) / sigma0(k);% nFeatures by 1
            theta = mu + noise * chol(Sigma);% nFeatures by 1
        end
        % Obtain a function sampled from the posterior GP.
        obj(:,k) = (theta' * sqrt(2 * sigma(k) / nFeatures) * cos(W * x' + repmat(b, 1, size(x, 1))))';
end
Cons = [];
end

%% MOMES calculation
function obj = evaluateMOMES(x, PFs, xx, yy, KernelMatrixInv, l, sigma, sigma0)
nx = length(x(:,1));
[nM,nK] = size(PFs);% equal to number of approximations
N_obj = length(yy(1,:));
objtt = zeros(nx,nM);
for i=1:nM
    u = zeros(nx,N_obj);  mse = zeros(nx,N_obj);
    for k = 1 : N_obj
        [u(:, k),mse(:, k)] = mean_var(x, xx, yy(:,k), KernelMatrixInv{k,i}, l{k}(i,:), sigma{k}(i), sigma0{k}(i));
        mse(mse(:, k)<=sigma0{k}(i)+eps, k) = sigma0{k}(i)+eps;
    end
    smse = sqrt(mse);  
    ysj_star = zeros(nK,N_obj);    objt = zeros(nx,nK);
    for s=1:nK
        ysj_star(s,:) = max(PFs{i,s});
        gamma_sj = (repmat(ysj_star(s,:),nx,1) - u)./smse;
        objt(:,s) = sum(gamma_sj.*Gaussian_PDF(gamma_sj)/2./Gaussian_CDF(gamma_sj)  ...
            - log(Gaussian_CDF(gamma_sj)), 2);
    end
    % the objective is to be maximized ???
    objtt(:,i) = sum(objt,2)/nK;
end
% the MOMES is to be maximized
obj = sum(objtt,2)/nM;
end

%% standard PDF of Gaussian
function y=Gaussian_PDF(x)
    y=1/sqrt(2*pi)*exp(-x.^2/2);
end

%% standard CDF of Gaussian
function y=Gaussian_CDF(x)
	y=0.5*(1+erf(x/sqrt(2)));
end