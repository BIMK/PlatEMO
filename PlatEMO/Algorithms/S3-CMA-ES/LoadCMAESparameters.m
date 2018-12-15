function tempPara = LoadCMAESparameters(Global, groups, popSize)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Huangke Chen

   para.lambda = popSize;	% population size, offspring number   
   tempPara    = repmat(para,1,length(groups));
   
   for g = 1 : length(groups)	% Set the parameters for each variable group
       
       DVindex = groups{g};     % the variable index of this group
       N = length(DVindex);     % variable number
       
       % set the CMA-ES parameters for each group
       tempPara(g).N = N;
       
       mu = tempPara(g).lambda/2;           % number of parents/points for recombination
       weights = log(mu+1/2) - log(1:mu)';  % muXone array for weighted recombination
       tempPara(g).mu = floor(mu);
       
       tempPara(g).weights = weights/sum(weights);	% normalize recombination weights array
       mueff = sum(weights)^2/sum(weights.^2);      % variance-effectiveness of sum w_i x_i
       tempPara(g).mueff = mueff;
       
       tempPara(g).cc    = (4+mueff/N)/(N+4+2*mueff/N);     % time constant for cumulation for C
       tempPara(g).cs    = (mueff+2)/(N+mueff+5);           % t-const for cumulation for sigma control
       tempPara(g).c1    = 2/((N+1.3)^2+mueff);             % learning rate for rank-one update of C
       tempPara(g).cmu   = min(1-tempPara(g).c1,2*(mueff-2+1/mueff)/((N+2)^2+mueff));	% and for rank-mu update
       tempPara(g).damps = 1 + 2*max(0,sqrt((mueff-1)/(N+1))-1) + tempPara(g).cs;       % damping for sigma
       tempPara(g).chiN  = N^0.5*(1-1/(4*N)+1/(21*N^2));	% expectation of ||N(0,I)|| == norm(randn(N,1))
       
       % variable parameters
       tempPara(g).xmean = Global.lower(DVindex)' + (Global.upper(DVindex)'-Global.lower(DVindex)').*rand(N,1);	% 10*rand(N, 1); objective variables initial point
       tempPara(g).sigma = 0.5;             % coordinate wise standard deviation (step size)
       tempPara(g).pc    = zeros(N,1);      % evolution paths for C and sigma
       tempPara(g).ps    = zeros(N,1);      % evolution paths for C and sigma
       B = eye(N,N);                        % B defines the coordinate system
       tempPara(g).B = B;
       D = ones(N,1);                       % diagonal D defines the scaling
       tempPara(g).D = D;
       tempPara(g).C = B*diag(D.^2)*B';     % covariance matrix C
       tempPara(g).eigeneval = 0;
       tempPara(g).counteval = 0;
   end
end