function [para,pop,BestVal,BestIndividual] = CMAES(para,bestmem,dim_index)  
% CMA-ES

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Huangke Chen

    N     = para.N;     % number of decision variables/problem dimension
    xmean = para.xmean;	% decision variables initial point
    sigma = para.sigma;	% coordinate wise standard deviation (step size)

    % Strategy parameter setting: Selection  
    lambda  = para.lambda;	% population size, offspring number
    mu      = para.mu;   	% floor(mu);        
    weights = para.weights;	% normalization of recombination weights array
    mueff   = para.mueff;	% variance-effectiveness of sum w_i x_i

    % Strategy parameter setting: Adaptation
    cc    = para.cc;   	% time constant for cumulation for C
    cs    = para.cs;  	% t-const for cumulation for sigma control
    c1    = para.c1;  	% learning rate for rank-one update of C
    cmu   = para.cmu;
    damps = para.damps;	% damping for sigma 
    chiN  = para.chiN;	% expectation of ||N(0,I)|| == norm(randn(N,1))

    % Initialize dynamic (internal) strategy parameters and constants
    pc        = para.pc;                % zeros(N, 1); 
    ps        = para.ps;                % zeros(N, 1); evolution paths for C and sigma
    B         = para.B;                 % eye(N, N); B defines the coordinate system
    D         = para.D;                 % ones(N, 1); diagonal D defines the scaling
    C         = para.C;                 % B * diag(D.^2) * B'; covariance matrix C
    invsqrtC  = B * diag(D.^-1) * B';	% C^-1/2 
    eigeneval = para.eigeneval;         % 0; track update of B and D  
    counteval = para.counteval;

    arx         = repmat(xmean,1,lambda) + sigma*B*(repmat(D,1,lambda).*randn(N,lambda));
    arx(arx<0)  = 0;
    arx(arx>10) = 10;

    tempDecs   = repmat(bestmem,lambda,1);
    pop        = arx';
    tempDecs(:,dim_index) = pop;
    population = INDIVIDUAL(tempDecs);
    objs       = population.objs;

    arfitness = sum(objs.^2,2);
    counteval = counteval + lambda;

    % Sort by fitness and compute weighted mean into xmean
    [arfitness,arindex] = sort(arfitness);	% minimization
    xold  = xmean;
    xmean = arx(:,arindex(1:mu))*weights;	% recombination, new mean value

    % Cumulation: Update evolution paths
    ps   = (1-cs)*ps ...
           + sqrt(cs*(2-cs)*mueff)*invsqrtC*(xmean-xold)/sigma;	% Eq.£¨24£©
    hsig = sum(ps.^2)/(1-(1-cs)^(2*counteval/lambda))/N < 2+4/(N+1);
    pc   = (1-cc)*pc ...
           + hsig*sqrt(cc*(2-cc)*mueff)*(xmean-xold)/sigma;

    % Adapt covariance matrix C
    artmp = (1/sigma)*(arx(:,arindex(1:mu))-repmat(xold,1,mu));	% mu difference vectors
    C     = (1-c1-cmu)*C ...                 	% regard old matrix
            + c1*(pc*pc' ...                 	% plus rank one update
            +(1-hsig)*cc*(2-cc)*C) ...      	% minor correction if hsig==0
            + cmu*artmp*diag(weights)*artmp';	% Eq.£¨30£© plus rank mu update

    % Adapt step size sigma
    sigma = sigma*exp((cs/damps)*(norm(ps)/chiN-1));

    % Update B and D from C
    if counteval - eigeneval > lambda/(c1+cmu)/N/10	% to achieve O(N^2)
        eigeneval = counteval;
        C = triu(C) + triu(C,1)';	% enforce symmetry
        if any(any(isnan(C))) || any(any(C==Inf))
            C = B*diag(D.^2)*B';
            C = triu(C) + triu(C,1)';	% enforce symmetry 
        end
        [B,D] = eig(C);	% eigen decomposition, B==normalized eigenvectors
        D     = sqrt(max(diag(D),0));     
    end
    BestVal        = arfitness(1);
    BestIndividual = population(arindex(1));

    % record the variable parameters
    para.xmean = xmean;     % objective variables initial point
    para.sigma = sigma;     % coordinate wise standard deviation (step size)
    para.pc    = pc;        % evolution paths for C and sigma
    para.ps    = ps;        % evolution paths for C and sigma
    para.B     = B;
    para.D     = D;
    para.C     = C;         % B * diag(D.^2) * B'; covariance matrix C
    para.eigeneval = eigeneval;
    para.counteval = counteval;
end