function [e, edata, eprior] = gp_e(w, gp, x, y, varargin)
%GP_E  Evaluate the energy function (un-normalized negative 
%      log marginal posterior)
%
%  Description
%    E = GP_E(W, GP, X, Y, OPTIONS) takes a Gaussian process
%    structure GP together with a matrix X of input vectors and a
%    matrix Y of targets, and evaluates the energy function E. Each
%    row of X corresponds to one input vector and each row of Y
%    corresponds to one target vector.
%
%    [E, EDATA, EPRIOR] = GP_E(W, GP, X, Y, OPTIONS) also returns
%    the data and prior components of the total energy. EDATA is
%    the negative marginal likelihood of the model.
%
%    The energy is minus log posterior cost function:
%        E = EDATA + EPRIOR 
%          = - log p(Y|X, th) - log p(th),
%    where th represents the parameters (lengthScale,
%    magnSigma2...), X is inputs and Y is observations (regression)
%    or latent values (non-Gaussian likelihood).
%
%    OPTIONS is optional parameter-value pair
%      z - optional observed quantity in triplet (x_i,y_i,z_i)
%          Some likelihoods may use this. For example, in case of
%          Poisson likelihood we have z_i=E_i, that is, expected
%          value for ith case.
%
%  See also
%    GP_G, GPCF_*, GP_SET, GP_PAK, GP_UNPAK
%
% Copyright (c) 2006-2010 Jarno Vanhatalo
% Copyright (c) 2010-2011 Aki Vehtari
% Copyright (c) 2010 Heikki Peura
% Copyright (c) 2014 Arno Solin and Jukka Koskenranta

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

if ~all(isfinite(w(:)));
  % instead of stopping to error, return NaN
  e=NaN;
  edata = NaN;
  eprior = NaN;
  return;
end

if isfield(gp,'latent_method') && ~strcmp(gp.latent_method,'MCMC')
  % use an inference specific method
  fh_e = gp.fh.e;
  switch nargout 
    case {0 1}
      [e] = fh_e(w, gp, x, y, varargin{:});
    case 2
      [e, edata] = fh_e(w, gp, x, y, varargin{:});
    case 3
      [e, edata, eprior] = fh_e(w, gp, x, y, varargin{:});
  end
  if ~isreal(e)
    warning('Energy is not real')
    e=NaN;
  end
  return
end

ip=inputParser;
ip.FunctionName = 'GP_E';
ip.addRequired('w', @(x) isempty(x) || isvector(x) && isreal(x));
ip.addRequired('gp',@isstruct);
ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addParamValue('z', [], @(x) isreal(x) && all(isfinite(x(:))))
ip.parse(w, gp, x, y, varargin{:});
z=ip.Results.z;

gp=gp_unpak(gp, w);
ncf = length(gp.cf);
n=size(x,1);
multicf=false;

if isfield(gp.lik, 'nondiagW') 
  % Help parameters for likelihoods with non-diagonal Hessian
  switch gp.lik.type
    case {'LGP', 'LGPC'}
      % Do nothing
    case {'Softmax', 'Multinom'}
      if size(y,1)~=size(x,1)
        y=reshape(y,size(x,1),size(y,1)/size(x,1));
      end
      [n,nout]=size(y);
      nl=[0 repmat(n,1,nout)];
      nl=cumsum(nl);
    otherwise
      if ~isfield(gp.lik,'xtime') && size(y,1)~=size(x,1)
        y=reshape(y,size(x,1),size(y,1)/size(x,1));
      end
      n=size(y,1);
      nout=length(gp.comp_cf);
      
      % Indices for looping over latent processes
      if ~isfield(gp.lik, 'xtime')
        nl=[0 repmat(n,1,nout)];
        nl=cumsum(nl);
      else
        xtime=gp.lik.xtime;
        ntime = size(xtime,1);
        n=n-ntime;
        nl=[0 ntime n];
        nl=cumsum(nl);
      end
  end
  if isfield(gp, 'comp_cf')  % own covariance for each ouput component
    multicf = true;
    if length(gp.comp_cf) ~= nout
      error('GP2_E: the number of component vectors in gp.comp_cf must be the same as number of outputs.')
    end
  end
end

% First Evaluate the data contribution to the error
switch gp.type
  % ============================================================
  % FULL GP (and compact support GP)
  % ============================================================
  case 'FULL'   % A full GP
    [K, C] = gp_trcov(gp, x);
    
    % Are there specified mean functions
    if  ~isfield(gp,'meanf')       % a zero mean function
      if issparse(C)            % compact support covariances are in use
        [LD,notpositivedefinite] = ldlchol(C);
        if notpositivedefinite
          [edata, eprior, e] = set_output_for_notpositivedefinite;
          return
        end
        edata = 0.5*(n.*log(2*pi) + sum(log(diag(LD))) + y'*ldlsolve(LD,y));
      else
        if ~isfield(gp.lik, 'nondiagW') || ismember(gp.lik.type, {'LGP' 'LGPC'})
          [L,notpositivedefinite] = chol(C,'lower');
          if notpositivedefinite
            [edata, eprior, e] = set_output_for_notpositivedefinite;
            return
          end
          ws=warning('off','MATLAB:singularMatrix');
          b=L\y;
          warning(ws);
          zc=sum(log(diag(L)));
        else
          b=zeros(nl(end),1);
          zc=0;
          y=y(:);
          if multicf
            for i1=1:nout
              if i1==1 && isfield(gp.lik, 'xtime')
                [tmp, C] = gp_trcov(gp, xtime, gp.comp_cf{i1});
              else
                [tmp, C] = gp_trcov(gp, x, gp.comp_cf{i1});
              end
              [L,notpositivedefinite]=chol(C,'lower');
              if notpositivedefinite
                [e, edata, eprior] = set_output_for_notpositivedefinite();
                return
              end
              %         b(:,i1) = L\y(:,i1);
              b(nl(i1)+1:nl(i1+1)) = L\y(nl(i1)+1:nl(i1+1));
              zc = zc + sum(log(diag(L)));
            end
          else
            [tmp, C] = gp_trcov(gp, x);
            [L,notpositivedefinite]=chol(C,'lower');
            if notpositivedefinite
              [e, edata, eprior] = set_output_for_notpositivedefinite();
              return
            end
            for i1=1:nout
              b(nl(i1)+1:nl(i1+1)) = L\y(nl(i1)+1:nl(i1+1));
              zc = zc + sum(log(diag(L)));
            end
          end
        end
        edata = 0.5*n.*log(2*pi) + zc + 0.5*b'*b;
      end
    else
      [H,b,B]=mean_prep(gp,x,[]);
      if isempty(C)
        L=1;
        C=0;
        logK=0;
        KH=H';
      elseif issparse(C)  
        [L,notpositivedefinite] = ldlchol(C);
        if notpositivedefinite
          [edata, eprior, e] = set_output_for_notpositivedefinite;
          return
        end
        logK = 0.5*sum(log(diag(L)));
        KH = ldlsolve(L,H');
      else          
        [L,notpositivedefinite] = chol(C,'lower');
        if notpositivedefinite
          [edata, eprior, e] = set_output_for_notpositivedefinite;
          return
        end
        logK = sum(log(diag(L)));
        KH = L'\(L\H');
      end
      
      A = B\eye(size(B)) + H*KH;
      [LA, notpositivedefinite] = chol(A,'lower');
      if notpositivedefinite
        [edata, eprior, e] = set_output_for_notpositivedefinite;
        return
      end
      M = H'*b-y;
      % When using CS-covariance function use matrix inversion lemma to
      % calculate the inverse -> faster computation
      if issparse(C)
          MNM = LA\(KH'*M);
          MNM = M'*ldlsolve(L,M) - MNM'*MNM;
      else
          [LN, notpositivedefinite] = chol(C + H'*B*H,'lower');
          if notpositivedefinite
              [edata, eprior, e] = set_output_for_notpositivedefinite;
              return
          end
          MNM = LN\M;
          MNM = MNM'*MNM;
      end
      
      edata = 0.5*MNM + logK + sum(log(diag(chol(B)))) + sum(log(diag(LA))) + 0.5*n*log(2*pi);
      
%       A = B\eye(size(B)) + H*KH;
%       M = H'*b-y;
%       [LN, notpositivedefinite] = chol(C + H'*B*H,'lower');
%       if notpositivedefinite
%         [edata, eprior, e] = set_output_for_notpositivedefinite;
%         return
%       end
%       MNM = LN\M;
%       MNM = MNM'*MNM;
%       
%       edata = 0.5*MNM + logK + 0.5*log(det(B)) + 0.5*log(det(A)) + 0.5*n*log(2*pi);
      
    end
    
    % ============================================================
    % FIC
    % ============================================================
  case 'FIC'
    % The equations in FIC are implemented as by Lawrence (2007)
    % See also Snelson and Ghahramani (2006) and Vanhatalo and Vehtari (2007)

    % First evaluate needed covariance matrices
    % v defines that parameter is a vector
    u = gp.X_u;
    [Kv_ff, Cv_ff] = gp_trvar(gp, x);  % n x 1  vector
    K_fu = gp_cov(gp, x, u);         % n x m
    K_uu = gp_trcov(gp, u);          % m x m, noiseles covariance K_uu
    K_uu = (K_uu+K_uu')./2;          % ensure the symmetry of K_uu
    [Luu, notpositivedefinite] = chol(K_uu,'lower');
    if notpositivedefinite
      [edata, eprior, e] = set_output_for_notpositivedefinite;
      return
    end
    % Evaluate the Lambda (La)
    % Q_ff = K_fu*inv(K_uu)*K_fu'
    % Here we need only the diag(Q_ff), which is evaluated below
    B=Luu\(K_fu');       % m x n
    Qv_ff=sum(B.^2)';
    Lav = Cv_ff-Qv_ff;   % n x 1, Vector of diagonal elements
                         % iLaKfu = diag(iLav)*K_fu = inv(La)*K_fu
    iLaKfu = zeros(size(K_fu));  % f x u,
    for i=1:n
      iLaKfu(i,:) = K_fu(i,:)./Lav(i);  % f x u
    end
    % The data contribution to the error is
    % E = n/2*log(2*pi) + 0.5*log(det(Q_ff+La)) + 0.5*y'inv(Q_ff+La)*y
    %   = + 0.5*log(det(La)) + 0.5*trace(iLa*y*y') - 0.5*log(det(K_uu))
    %     + 0.5*log(det(A)) - 0.5*trace(inv(A)*iLaKfu'*y*y'*iLaKfu)

    % First some help matrices...
    % A = chol(K_uu+K_uf*inv(La)*K_fu))
    A = K_uu+K_fu'*iLaKfu;
    A = (A+A')./2;     % Ensure symmetry
    [A, notpositivedefinite] = chol(A,'upper');
    if notpositivedefinite
      [edata, eprior, e] = set_output_for_notpositivedefinite;
      return
    end
    % The actual error evaluation
    % 0.5*log(det(K)) = sum(log(diag(L))), where L = chol(K). NOTE! chol(K) is upper triangular
    b = (y'*iLaKfu)/A;
    edata = sum(log(Lav)) + y'./Lav'*y - 2*sum(log(diag(Luu))) + 2*sum(log(diag(A))) - b*b';
    edata = .5*(edata + n*log(2*pi));
    % ============================================================
    % PIC
    % ============================================================
  case {'PIC' 'PIC_BLOCK'}
    % First evaluate needed covariance matrices
    % v defines that parameter is a vector
    u = gp.X_u;
    ind = gp.tr_index;
    K_fu = gp_cov(gp, x, u);         % n x m
    K_uu = gp_trcov(gp, u);    % m x m, noiseles covariance K_uu
    K_uu = (K_uu+K_uu')./2;     % ensure the symmetry of K_uu
    [Luu, notpositivedefinite] = chol(K_uu,'lower');
    if notpositivedefinite
      [edata, eprior, e] = set_output_for_notpositivedefinite;
      return
    end

    % Evaluate the Lambda (La)
    % Q_ff = K_fu*inv(K_uu)*K_fu'
    % Here we need only the blockdiag(Q_ff), which is evaluated below
    B=Luu\(K_fu');       % u x f  and B'*B = K_fu*K_uu*K_uf
    iLaKfu = zeros(size(K_fu));  % f x u
    edata = 0;
    for i=1:length(ind)        
      Qbl_ff = B(:,ind{i})'*B(:,ind{i});
      [Kbl_ff, Cbl_ff] = gp_trcov(gp, x(ind{i},:));
      Labl{i} = Cbl_ff - Qbl_ff;
      iLaKfu(ind{i},:) = Labl{i}\K_fu(ind{i},:);
      [Ltmp, notpositivedefinite]=chol(Labl{i},'upper');
      if notpositivedefinite
        [edata, eprior, e] = set_output_for_notpositivedefinite;
        return
      end
      edata = edata + 2*sum(log(diag(Ltmp))) + y(ind{i},:)'*(Labl{i}\y(ind{i},:));
    end
    % The data contribution to the error is
    % E = n/2*log(2*pi) + 0.5*log(det(Q_ff+La)) + 0.5*y'inv(Q_ff+La)y

    % First some help matrices...
    % A = chol(K_uu+K_uf*inv(La)*K_fu))
    A = K_uu+K_fu'*iLaKfu;
    A = (A+A')./2;     % Ensure symmetry
    [A, notpositivedefinite] = chol(A,'lower');
    if notpositivedefinite
      [edata, eprior, e] = set_output_for_notpositivedefinite;
      return
    end
    % The actual error evaluation
    % 0.5*log(det(K)) = sum(log(diag(L))), where L = chol(K). NOTE! chol(K) is upper triangular
    b = (y'*iLaKfu)*inv(A)';
    edata = edata - 2*sum(log(diag(Luu))) + 2*sum(log(diag(A))) - b*b';
    edata = .5*(edata + n*log(2*pi));
    % ============================================================
    % CS+FIC
    % ============================================================
  case 'CS+FIC'
    u = gp.X_u;

    % Separate the FIC and CS covariance functions
    cf_orig = gp.cf;

    cf1 = {};
    cf2 = {};
    j = 1;
    k = 1;
    for i = 1:ncf
      if ~isfield(gp.cf{i},'cs')
        cf1{j} = gp.cf{i};
        j = j + 1;
      else
        cf2{k} = gp.cf{i};
        k = k + 1;
      end
    end
    gp.cf = cf1;

    % Evaluate the covariance matrices needed for FIC part
    [Kv_ff, Cv_ff] = gp_trvar(gp, x);  % n x 1  vector
    K_fu = gp_cov(gp, x, u);         % n x m
    K_uu = gp_trcov(gp, u);    % m x m, noiseles covariance K_uu
    K_uu = (K_uu+K_uu')./2;     % ensure the symmetry of K_uu
    [Luu, notpositivedefinite] = chol(K_uu,'lower');
    if notpositivedefinite
      [edata, eprior, e] = set_output_for_notpositivedefinite;
      return
    end

    % Evaluate the Lambda (La)
    % Q_ff = K_fu*inv(K_uu)*K_fu'
    B=Luu\(K_fu');       % m x n
    Qv_ff=sum(B.^2)';
    Lav = Cv_ff-Qv_ff;   % n x 1, Vector of diagonal elements

    % Evaluate the CS covariance matrix
    gp.cf = cf2;
    K_cs = gp_trcov(gp,x);
    La = sparse(1:n,1:n,Lav,n,n) + K_cs;
    gp.cf = cf_orig;     % Set the original covariance functions in the GP structure

    [LD, notpositivedefinite] = ldlchol(La);
    if notpositivedefinite
      [edata, eprior, e] = set_output_for_notpositivedefinite;
      return
    end
    
    %        iLaKfu = La\K_fu;
    iLaKfu = ldlsolve(LD,K_fu);
    edata = sum(log(diag(LD))) + y'*ldlsolve(LD,y);
    % The data contribution to the error is
    % E = n/2*log(2*pi) + 0.5*log(det(Q_ff+La)) + 0.5*y'inv(Q_ff+La)y

    % First some help matrices...
    % A = chol(K_uu+K_uf*inv(La)*K_fu))
    A = K_uu+K_fu'*iLaKfu;
    A = (A+A')./2;     % Ensure symmetry
    [A, notpositivedefinite] = chol(A,'upper');
    if notpositivedefinite
      [edata, eprior, e] = set_output_for_notpositivedefinite;
      return
    end
    % The actual error evaluation
    % 0.5*log(det(K)) = sum(log(diag(L))), where L = chol(K). NOTE! chol(K) is upper triangular
    %b = (y'*iLaKfu)*inv(A)';
    b = (y'*iLaKfu)/A;
    edata = edata - 2*sum(log(diag(Luu))) + 2*sum(log(diag(A))) - b*b';
    edata = .5*(edata + n*log(2*pi));
    % ============================================================
    % DTC/VAR
    % ============================================================
  case {'DTC' 'VAR' 'SOR'}
    % Implementation of DTC varies only slightly from FIC: essentially, only
    % Lav is defined differently. For equations, see e.g. Quinonero-Candela
    % and Rasmussen. For VAR, a trace term is added to the DTC model, see 
    % Titsias (2009).
    
    % First evaluate needed covariance matrices
    % v defines that parameter is a vector
    u = gp.X_u;
    [Kv_ff, Cv_ff] = gp_trvar(gp, x);  % n x 1  vector
    K_fu = gp_cov(gp, x, u);         % n x m
    K_uu = gp_trcov(gp, u);          % m x m, noiseles covariance K_uu
    K_uu = (K_uu+K_uu')./2;          % ensure the symmetry of K_uu
    [Luu, notpositivedefinite] = chol(K_uu, 'lower');
    if notpositivedefinite
      [edata, eprior, e] = set_output_for_notpositivedefinite;
      return
    end
    % Evaluate the Lambda (La)
    % Q_ff = K_fu*inv(K_uu)*K_fu';
    % Here we need only the diag(Q_ff), which is evaluated below
    B=Luu\(K_fu');       % m x n
    Qv_ff=sum(B.^2)';
    Lav = Cv_ff-Kv_ff;   % n x 1, Vector of diagonal elements
                         % iLaKfu = diag(iLav)*K_fu = inv(La)*K_fu
    iLaKfu = zeros(size(K_fu));  % f x u,
    for i=1:n
      iLaKfu(i,:) = K_fu(i,:)./Lav(i);  % f x u
    end
    % The data contribution to the error is
    % E = n/2*log(2*pi) + 0.5*log(det(Q_ff+La)) + 0.5*t'inv(Q_ff+La)*t
    %   = + 0.5*log(det(La)) + 0.5*trace(iLa*t*t') - 0.5*log(det(K_uu))
    %     + 0.5*log(det(A)) - 0.5*trace(inv(A)*iLaKfu'*t*t'*iLaKfu)

    % First some help matrices...
    % A = chol(K_uu+K_uf*inv(La)*K_fu))
    A = K_uu+K_fu'*iLaKfu;
    A = (A+A')./2;     % Ensure symmetry
    [A, notpositivedefinite] = chol(A);
    if notpositivedefinite
      [edata, eprior, e] = set_output_for_notpositivedefinite;
      return
    end
    % The actual error evaluation
    % 0.5*log(det(K)) = sum(log(diag(L))), where L = chol(K). NOTE! chol(K) is upper triangular
    b = (y'*iLaKfu)/A;
    edata = sum(log(Lav)) + y'./Lav'*y - 2*sum(log(diag(Luu))) + 2*sum(log(diag(A))) - b*b';
    edata = 0.5*(edata + n*log(2*pi));
    if strcmp(gp.type, 'VAR')
      edata = edata + 0.5*sum((Kv_ff-Qv_ff)./Lav);
    end
    %edata = edata - 0.5*sum((Kv_ff-Qv_ff)./Lav);% - sum(diag(B'*B),1)); %sum(B.^2,1)'
    %sum(Qv_ff)
    %K_ff=gp_trcov(gp,x);
    %0.5*trace(K_ff-K_fu*inv(K_uu)*K_fu')
    %0.5*trace(K_ff-B'*B)
    
  case {'KALMAN'}  
    % ============================================================
    % Kalman filtering and smoothing
    % ============================================================
    %
    % The implementation below is primarily based on the methods 
    % presented in the following publication. If you find this 
    % useful as a part of your own research, please cite the papers.
    %
    %  [1] Simo Sarkka, Arno Solin, Jouni Hartikainen (2013).
    %      Spatiotemporal learning via infinite-dimensional Bayesian
    %      filtering and smoothing. IEEE Signal Processing Magazine,
    %      30(4):51-61.
    %
    %  [2] Simo Sarkka (2013). Bayesian filtering and smoothing. 
    %      Cambridge University Press.
    %
    %  [3] Simo Sarkka (2006). Recursive Bayesian inference on stochastic
    %      differential equations. Doctoral dissertation, Helsinki 
    %      University of Technology, Filand.
    %
    
    % Ensure that this is a purely temporal problem
    if size(x,2) > 1,
      error('The ''KALMAN'' option only supports one-dimensional data.')  
    end
    
    % Extract the noise magnitude from the GP likelihood model
    R = gp.lik.sigma2;
    
    % Initialize model matrices
    F    = [];
    L    = [];
    Qc   = [];
    H    = [];
    Pinf = [];
    
    % For each covariance function
    for j=1:length(gp.cf)
 
      % Form state-space model from the gp.cf{j}
      try
        [jF,jL,jQc,jH,jPinf] = gp.cf{j}.fh.cf2ss(gp.cf{j});
      catch
        [edata, eprior, e] = set_output_for_notpositivedefinite;
        return 
      end
        
      % Stack model
      F    = blkdiag(F,jF);
      L    = blkdiag(L,jL);
      Qc   = blkdiag(Qc,jQc);
      H    = [H jH];    
      Pinf = blkdiag(Pinf,jPinf);
      
    end
    
    % Sort values
    [x,ind] = sort(x);
    y = y(ind);
    
    % Set initial state
    m = zeros(size(F,1),1);
    P = Pinf;

    % Initialize likelihood
    edata = 0;
    
    % Run filter for evaluating the marginal likelihood
    for k=1:size(y,1)
        
      % Solve A using the method by Davison
      if (k>1)
            
        % Discrete-time solution (only for stable systems)
        dt = x(k)-x(k-1);
        A  = expm(F*dt);
        Q  = Pinf - A*Pinf*A';
        
        % Prediction step
        m = A * m;
        P = A * P * A' + Q;
        
      end
      
      % Update step
      S = H*P*H'+R;
      K = P*H'/S;
      v = y(k)-H*m;
      m = m + K*v;
      P = P - K*H*P;
      
      % Update log likelihood
      edata = edata + 1/2*log(det(2*pi*S)) + 1/2*v'/S*v;
      
    end
    
  otherwise
    error('Unknown type of Gaussian process!')
end

% ============================================================
% Evaluate the prior contribution to the error from covariance functions
% ============================================================
eprior = 0;
if ~isempty(strfind(gp.infer_params, 'covariance'))
  for i=1:ncf
    gpcf = gp.cf{i};
    eprior = eprior - gpcf.fh.lp(gpcf);
  end
end

% ============================================================
% Evaluate the prior contribution to the error from Gaussian likelihood
% ============================================================
if ~isempty(strfind(gp.infer_params, 'likelihood')) && isfield(gp.lik.fh,'trcov') && isfield(gp.lik, 'p')
  % a Gaussian likelihood
  lik = gp.lik;
  eprior = eprior - lik.fh.lp(lik);
end

% ============================================================
% Evaluate the prior contribution to the error from the inducing inputs
% ============================================================
if ~isempty(strfind(gp.infer_params, 'inducing'))
  if isfield(gp, 'p') && isfield(gp.p, 'X_u') && ~isempty(gp.p.X_u)
    if iscell(gp.p.X_u) % Own prior for each inducing input
      for i = 1:size(gp.X_u,1)        
        pr = gp.p.X_u{i};
        eprior = eprior - pr.fh.lp(gp.X_u(i,:), pr);        
      end
    else
      eprior = eprior - gp.p.X_u.fh.lp(gp.X_u(:), gp.p.X_u);
    end
  end
end

% ============================================================
% Evaluate the prior contribution to the error from mean functions
% ============================================================
if ~isempty(strfind(gp.infer_params, 'mean'))
  for i=1:length(gp.meanf)
    gpmf = gp.meanf{i};
    eprior = eprior - gpmf.fh.lp(gpmf);
  end
end

e = edata + eprior;

function [edata, eprior, e] = set_output_for_notpositivedefinite()
  %instead of stopping to chol error, return NaN
  edata = NaN;
  eprior= NaN;
  e = NaN;
end

end

