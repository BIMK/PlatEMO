function [P_temp, P_temp_rank] = OSP(ger_end, Pt, Pt_apt, Pt_Rank, Problem, ger_init, p)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Elaine Guerrero-Pena

    NP = Problem.N;
    m = Problem.M;

    intv_ger = ger_end - ger_init;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data pre-processing

    phi_data = [];

    % Current generation data
    C = ger_end;
    S = find(Pt_Rank(:,C)==1);
    phi_apt = Pt_apt(S,:,C);
    indiv_M = Pt(S,:,C);
    sol_phiS = sum(Pt_Rank(:,ger_end)==1);

    [~, idx_sortmax] = sort(phi_apt(:,1));
    phi_S = phi_apt(idx_sortmax,:);
    indiv_S = indiv_M(idx_sortmax,:);

    for j = 1:intv_ger+1
        idx = []; idx1 = [];

        ger_i = j + ger_init-1;
        idx1 = find(Pt_Rank(:,ger_i)==1);

        if length(idx1) ~= sol_phiS
            [~,idx] = min(pdist2(phi_S,Pt_apt(idx1,:,ger_i)),[],2);
        else
            [~, idx] = sort(Pt_apt(idx1,1,ger_i));
        end

        phi_data(:,:,j) = Pt_apt(idx1(idx),:,ger_i);
    end

    % Find the "best individual"
    Dist_apt = pdist2(phi_data(:,:,1),phi_data(:,:,end));
    dist_apt = Dist_apt(sub2ind(size(Dist_apt),1:size(Dist_apt,1),1:size(Dist_apt,2)));
    idx_dif = find(sum((phi_data(:,:,end) - phi_data(:,:,1)) < 0,2)==m);
    [~,idx_apt] = max(dist_apt(idx_dif));
    best_indiv = idx_dif(idx_apt);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Forecasting of the front movement using ARX time series models
    phi_f = fun_ARX(phi_data,best_indiv,p);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Searching for the solutions in the decision space that have fitness values
    % in the objective space as close as possible to those of the predicted front
    X_f = fun_Otimiz(indiv_S,phi_S,phi_f,Problem);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Complementing a new population using Gaussian Mixture Model-based Local
    % Search (GMM-LS)
    ksol = NP - length(X_f);

    if ksol~=0
        X_LS = fun_GMMOV(X_f.decs, Problem, ksol);
        P_temp = [X_f X_LS];
    else
        P_temp = X_f;
    end

    P_temp_rank = NDSort(P_temp.objs,P_temp.cons,Problem.N);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Block 1: Forecast with ARX model
function F_p = fun_ARX(indiv_m,idx_model,p)

    warning off

    na = eye(size(indiv_m,2));
    nb = [];
    nk = [];

    if ~isempty(idx_model) % ARX model is built using the movement of best
        % individual in the generations
        y1 = reshape(indiv_m(idx_model,:,:),size(indiv_m,2),size(indiv_m,3))';
        y2 = y1(all(~isnan(y1),2),:);
        [~,idx_y] = unique(y2(:,1),'stable');
        y = y2(idx_y,:);
        if size(y,1) < 3
            y=y2;
        end

        data = iddata(y,[],'TimeUnit', 'hours');
        Mdl = arx(data, [na nb nk]); % Estimate ARX model parameters using the
        % best individual

        for i = 1:size(indiv_m,1)

            y1 = reshape(indiv_m(i,:,:),size(indiv_m,2),size(indiv_m,3))';
            y = y1(all(~isnan(y1),2),:);

            data = iddata(y,[],'TimeUnit', 'hours');

            ydata = forecast(Mdl,data,p); % Forecast with the ARX Model found

            F_p(i,:) = ydata.y(end,:);   % Predictive Front

        end

    else % ARX model is built for the individuals using the movement of each
        % solution in the generations
        for i = 1:size(indiv_m,1)

            y1 = reshape(indiv_m(i,:,:),size(indiv_m,2),size(indiv_m,3))';
            y = y1(all(~isnan(y1),2),:);

            data = iddata(y,[],'TimeUnit', 'hours');

            Mdl = arx(data, [na nb nk]); % Estimate ARX model parameters using
            % each solution

            ydata = forecast(Mdl,data,p); % Forecast with the ARX Model found

            F_p(i,:) = ydata.y(end,:);   % Predictive Front

        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Block 2: Leading solutions approximation
function XF = fun_Otimiz(indiv_Max,apt_Max,yF,Problem)

    for i=1:size(yF,1)

        x0 = indiv_Max(i,:);
        f_x0 = apt_Max(i,:);

        fun = @(x)pdist2(Problem.Evaluation(x).obj,yF(i,:));

        XF_temp  = Problem.Evaluation(x0);
        dist_apt = pdist2(f_x0,yF(i,:));

        % Intermedial optimization using sequential quadratic programming method (SQP)
        options = optimoptions('fmincon','Display','off','Algorithm','sqp');
        [x, dist_aptx] = fmincon(fun,x0,[],[],[],[],...
            Problem.lower,Problem.upper,[],options);

        if dist_aptx < dist_apt
            XF_temp = Problem.Evaluation(x);
        end

        XF(i) = XF_temp; % "leading individuals"
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Block 3: Complete the population using GMM-LS
function Y_LS = fun_GMMOV(indiv, Problem, ksol)

    lbound = Problem.lower;
    ubound = Problem.upper;

    % Determine GMM using Variational Inference
    GMModel = GMMModel_VI(indiv,size(indiv,1));

    % Sample ksol individuals following the GMM distribution
    Y = fun_GMMRand(GMModel,Problem.D,ksol);

    % Repair non-feasible solutions
    for iy = 1:size(Y,1)
        Y(iy,Y(iy,:) > ubound) =  ubound(Y(iy,:) > ubound);
        Y(iy,Y(iy,:) < lbound) =  lbound(Y(iy,:) < lbound);
    end

    % Evaluate new individuals
    Y_LS = Problem.Evaluation(Y);

end

% Generate samples from a GMM.
function [X, z] = fun_GMMRand(model,Dp,n)
    
    % GMM
    K = size(model.mean,2);
    w = model.weights;
    mu = model.mean;
    Sigma = model.Sigma;
    
    % Generate samples from a multinomial distribution.
    p = cumsum(w(:));
    [~,z] = histc(rand(1,n),[0;p/p(end)]);

    % Generate samples from a Gaussian distribution.
    X = zeros(n,Dp);
    
    for i = 1:K
        [R,err] = chol(Sigma(:,:,i));
        if err ~= 0
            R = Sigma(:,:,i);
        end
        X(z==i,:) = (R'*randn(size(R,1),sum(z==i))+repmat(mu(:,i),1,sum(z==i)))';
    end
    
end

%% Variational Inference
function GMModel = GMMModel_VI(indiv1,K_G)

    [n,d] = size(indiv1);

    prior.alpha = 0.1;
    prior.kappa = 0.001;
    prior.m = mean(indiv1)';
    prior.v = d;
    prior.M = 0.00001*eye(d);
    prior.logW = -2*sum(log(diag(chol(prior.M))));

    tol = 1e-8;
    MaxIter = 2000;
    LB = -inf(1,MaxIter);
    
    label = ceil(K_G*rand(1,n));
    GMModel.R = full(sparse(1:n,label,1,n,K_G,n));
    GMModel = fun_maximize(indiv1',GMModel,prior);

    for iter = 2:MaxIter
        GMModel = fun_expect(indiv1',GMModel);
        GMModel = fun_maximize(indiv1',GMModel,prior);
        LB(iter) = fun_bound(indiv1',GMModel,prior)/n;
        if abs(LB(iter)-LB(iter-1)) < tol*abs(LB(iter))
            break;
        end
    end

    for i = 1:K_G
        [~,errm] = chol(GMModel.Sigma(:,:,i));
        if errm ~= 0
            Sigmam = (GMModel.Sigma(:,:,i) + GMModel.Sigma(:,:,i).')/2;
            GMModel.Sigma(:,:,i) = Sigmam;
        end
    end

end

%% Maximize Step
function GMModel = fun_maximize(X, GMModel, prior)

    R = GMModel.R;
    m = bsxfun(@plus,prior.kappa*prior.m,X*R);

    k = size(m,2);
    r = sqrt(R');

    for i = 1:k
        X_m = bsxfun(@times,bsxfun(@minus,X,m(:,i)),r(i,:));
        p_m = prior.m-m(:,i);
        U(:,:,i) = chol(prior.M + X_m*X_m' + prior.kappa*(p_m*p_m'));
        log_W(i) = -2*sum(log(diag(U(:,:,i))));
    end

    GMModel.alpha = prior.alpha + sum(R,1);
    GMModel.kappa = prior.kappa + sum(R,1);
    GMModel.mean = bsxfun(@times,m,1./GMModel.kappa);
    GMModel.v = prior.v + sum(R,1);
    GMModel.Sigma = U;
    GMModel.logW = log_W;
end

%% Expectation Step
function GMModel = fun_expect(X, GMModel)

    [d,k] = size(GMModel.mean);

    for i = 1:k
        q = (GMModel.Sigma(:,:,i)'\bsxfun(@minus,X,GMModel.mean(:,i)));
        E_q(:,i) = d/GMModel.kappa(i) + GMModel.v(i)*dot(q,q,1);
    end
    E_logL = sum(psi(0,0.5*bsxfun(@minus,GMModel.v+1,(1:d)')),1) + ...
        d*log(2)+GMModel.logW;
    E_logpi = psi(0,GMModel.alpha) - psi(0,sum(GMModel.alpha));
    log_Rho = bsxfun(@plus,-0.5*bsxfun(@minus,E_q,E_logL - d*log(2*pi)),E_logpi);

    y = max(log_Rho,[],2);
    log_sumexp = y+log(sum(exp(bsxfun(@minus,log_Rho,y)),2));
    ii = isinf(y);
    if any(ii(:))
        log_sumexp(ii) = y(ii);
    end
    
    GMModel.logR = bsxfun(@minus,log_Rho,log_sumexp);
    R = exp(GMModel.logR);
    GMModel.R = R;

    GMModel.weights = zeros(1, size(GMModel.R,2));
    for i=1:size(GMModel.R, 2)
        GMModel.weights(1,i) = sum(GMModel.R(:,i)/size(GMModel.R,1));
    end

end

%% Bound
function LB = fun_bound(X, GMModel, prior)
d = size(X,1);
k = size(GMModel.R,2);

Ep_z = 0;
Eq_z = dot(GMModel.R(:),GMModel.logR(:));
Ep_pi = gammaln(k*prior.alpha)-k*gammaln(prior.alpha);
Eq_pi = gammaln(sum(GMModel.alpha))-sum(gammaln(GMModel.alpha));
Ep_mu = 0.5*d*k*log(prior.kappa);
Eq_mu = 0.5*d*sum(log(GMModel.kappa));
Ep_L = k*(-0.5*prior.v*(prior.logW+d*log(2))-fun_logMvGamma(0.5*prior.v,d));
Eq_L = sum(-0.5*GMModel.v.*(GMModel.logW+d*log(2))-fun_logMvGamma(0.5*GMModel.v,d));
Ep_X = -0.5*d*size(X,2)*log(2*pi);
LB = Ep_z-Eq_z+Ep_pi-Eq_pi+Ep_mu-Eq_mu+Ep_L-Eq_L+Ep_X;
end

function y = fun_logMvGamma(x,d)
dim = size(x);
x = reshape(x,1,prod(dim));
x = bsxfun(@plus,repmat(x,d,1),(1-(1:d)')/2);
y = d*(d-1)/4*log(pi)+sum(gammaln(x),1);
y = reshape(y,dim);
end