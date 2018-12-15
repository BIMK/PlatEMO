function [Population,Rank,zeta,zmin] = EnvironmentalSelection(Population,zeta,zmin,K)
% The environmental selection of one-by-one EA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    PopObj = Population.objs;
    [N,M]  = size(PopObj);
    Q      = 1 : N;         % Set of remaining solutions
    Qs     = [];            % Set of preselected solutions
    Qth    = [];            % Set of de-emphasized solutions
    Qd     = [];            % Set of dominated solutions
    Rank   = inf(1,N);      % Rank of each solution in one-by-one selection
    nrank  = 1;             % Current rank

    %% Normalization
    % The normalization is optional in this algorithm, thus no
    % normalization here due to the poor performance on DTLZ
    zmin   = min(zmin,min(PopObj,[],1));
    PopObj = PopObj - repmat(zmin,N,1);
    
    %% Convergence and distribution indicators of each solution
    c = sqrt(sum(PopObj.^2,2));
    d = pdist2(PopObj,PopObj,'cosine');

    %% Select the corner solutions
    for m = 1 : M
        if ~isempty(Q)
            [xunique,Q,Qs] = SelectingUniqueIndividual(sqrt(sum(PopObj(:,[1:m-1,m+1:end]).^2,2)),Q,Qs);
            [Q,Qth,Qd]     = DeemphasizingIndividual(xunique,zeta,Q,Qth,Qd,PopObj,d);
            Rank(xunique)  = nrank;
        end
    end
    
    %% One-by-one selection
    while ~isempty(Q)
        [xunique,Q,Qs] = SelectingUniqueIndividual(c,Q,Qs);
        [Q,Qth,Qd]     = DeemphasizingIndividual(xunique,zeta,Q,Qth,Qd,PopObj,d);
        Rank(xunique)  = nrank;
    end
    r  = length(Qs)/K;
    Nd = length(Qd);
    
    %% Make the number of selected solutions be K
    while length(Qs) < K
        if isempty(Q)
            Q     = [Qth,Qd];
            Qth   = [];
            Qd    = [];
            nrank = nrank + 1;
        end
        [xunique,Q,Qs] = SelectingUniqueIndividual(c,Q,Qs);
        [Q,Qth,Qd]     = DeemphasizingIndividual(xunique,zeta,Q,Qth,Qd,PopObj,d);
        Rank(xunique)  = nrank;
    end
    Population = Population(Qs(1:K));
	Rank       = Rank(Qs(1:K));
    
    %% Update the distribution threshold
    if Nd <= K
        zeta = zeta*exp((r-1)/M);
    end
end

function [xunique,Q,Qs] = SelectingUniqueIndividual(Fitness,Q,Qs)
% Select the solution having the least fitness value

    [~,x]   = min(Fitness(Q));
    xunique = Q(x);
    Qs      = [Qs,xunique];
    Q(x)    = [];
end

function [Q,Qth,Qd] = DeemphasizingIndividual(xunique,zeta,Q,Qth,Qd,PopObj,d)
% De-emphasize the neighboring or dominated solutions

    xi    = d(Q,xunique) < zeta;
    Qth   = [Qth,Q(xi)];
    Q(xi) = [];
    xi    = all(PopObj(Q,:)>=repmat(PopObj(xunique,:),length(Q),1),2);
    Qd    = [Qd,Q(xi)];
    Q(xi) = [];
end