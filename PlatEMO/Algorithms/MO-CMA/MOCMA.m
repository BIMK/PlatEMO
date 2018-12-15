function MOCMA(Global)
% <algorithm> <M>
% Multi-objective covariance matrix adaptation evolution strategy

%------------------------------- Reference --------------------------------
% C. Igel, N. Hansen, and S. Roth, Covariance matrix adaptation for multi-
% objective optimization, Evolutionary computation, 2007, 15(1): 1-28.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate the initial individuals in CMA-ES
    Population = Global.Initialization();
    ptarget    = 1/5.5;
    a          = struct('x',num2cell(Population.decs,2)','psucc',ptarget,'sigma',0.5,'pc',0,'C',eye(Global.D),'Individual',num2cell(Population));

    %% Optimization
    while Global.NotTermination([a.Individual])
        % Generate new individuals
        for k = 1 : Global.N
            a1(k)            = a(k);
            a1(k).x          = mvnrnd(a(k).x,a(k).sigma^2*a(k).C,1);
            a1(k).Individual = INDIVIDUAL(a1(k).x);
        end
        
        % Update the fitness of each individual
        Q           = [a,a1];
        Population  = [Q.Individual];
        % Penalized fitness for handling box constraints
        PopObj      = Population.objs + repmat(1e-6*sum((cat(1,Q.x)-Population.decs).^2,2),1,Global.M);
        % Calculate the fitness of each individual
        FrontNo     = NDSort(PopObj,inf);
        CrowdDis    = CrowdingDistance(PopObj,FrontNo);
        [~,rank]    = sortrows([FrontNo;-CrowdDis]');
        [~,fitness] = sort(rank);
        
        % Update the CMA models
        for k = 1 : Global.N
            a(k)  = updateStepSize(a(k),fitness(Global.N+k)<fitness(k),ptarget);
            a1(k) = updateStepSize(a1(k),fitness(Global.N+k)<fitness(k),ptarget);
            a1(k) = updateCovariance(a1(k),(a1(k).x-a(k).x)/a(k).sigma);
        end
        
        % Update Q
        Q = [a,a1];

        % Individuals for next generation
        Q = [a,a1];
        a = Q(rank(1:Global.N));
    end
end