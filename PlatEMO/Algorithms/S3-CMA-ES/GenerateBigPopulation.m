function [BigPopulation,tempPara] = GenerateBigPopulation(PV,groups,Archive)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Huangke Chen

   cDecs  = Archive.decs;
   PVDecs = cDecs(:,PV);
   
   BigPopulation = cell(1,size(PVDecs,1));	% initialize the big population to contain all the population
   cPopSize      = 16;                      % the population size of each population
   
   para.lambda = cPopSize;                                      % population size, offspring number
   tempPara    = repmat(para,size(PVDecs,1),length(groups));	% Initialize the parameters for CMA-ES
   
   for ci = 1: size(PVDecs,1)
       % Construct the population
       BigPopulation{ci} = repmat(cDecs(ci,:),cPopSize,1);
       
       % Set the CMA-ES parameters
       for g = 1 : length(groups)	% Set the parameters for each variable group
           DVindex = groups{g};     % the variable index of this group
           N = length(DVindex);     % variable number
           
           % set the CMA-ES parameters for each group
           tempPara(ci,g).N = N;
           mu      = cPopSize/2;                % number of parents/points for recombination
           weights = log(mu+1/2) - log(1:mu)';	% muXone array for weighted recombination
           tempPara(ci,g).mu = floor(mu);
           
           tempPara(ci,g).weights = weights/sum(weights);	% normalize recombination weights array
           mueff = sum(weights)^2/sum(weights.^2);          % variance-effectiveness of sum w_i x_i
           tempPara(ci,g).mueff = mueff;
           
           tempPara(ci,g).cc    = (4+mueff/N)/(N+4+2*mueff/N);                                  % time constant for cumulation for C
           tempPara(ci,g).cs    = (mueff+2)/(N+mueff+5);                                        % t-const for cumulation for sigma control
           tempPara(ci,g).c1    = 2/((N+1.3)^2+mueff);                                          % learning rate for rank-one update of C
           tempPara(ci,g).cmu   = min(1-tempPara(ci,g).c1,2*(mueff-2+1/mueff)/((N+2)^2+mueff));	% and for rank-mu update
           tempPara(ci,g).damps = 1 + 2*max(0,sqrt((mueff-1)/(N+1))-1) + tempPara(ci, g).cs;    % damping for sigma
           tempPara(ci,g).chiN  = N^0.5*(1-1/(4*N)+1/(21*N^2));                                 % expectation of ||N(0,I)|| == norm(randn(N,1))
           
           % variable parameters
           paraDecs = BigPopulation{ci};
           tempPara(ci,g).xmean = mean(paraDecs(:,DVindex))';
           tempPara(ci,g).sigma = 0.1;          % coordinate wise standard deviation (step size)
           tempPara(ci,g).pc    = zeros(N,1);	% evolution paths for C and sigma
           tempPara(ci,g).ps    = zeros(N,1);	% evolution paths for C and sigma
           B = eye(N,N);	% B defines the coordinate system
           tempPara(ci,g).B = B;
           D = ones(N,1);	% diagonal D defines the scaling
           tempPara(ci,g).D = D;
           tempPara(ci,g).C = B*diag(D.^2)*B';	% covariance matrix C
           tempPara(ci,g).eigeneval = 0;
           tempPara(ci,g).counteval = 0;
       end
   end
end