function MOEADCMA(Global)
% <algorithm> <M>
% MOEA/D with covariance matrix adaptation evolution strategy
% K --- 5 --- Number of groups

%------------------------------- Reference --------------------------------
% H. Li, Q. Zhang, and J. Deng, Biased multiobjective optimization and
% decomposition algorithm, IEEE Transactions on Cybernetics, 2017, 47(1):
% 52-66.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    K = Global.ParameterSet(5);

    %% Generate the weight vectors
    [W,Global.N] = UniformPoint(Global.N,Global.M);
    T = ceil(Global.N/10);
    % Transformation on W
    W = 1./W./repmat(sum(1./W,2),1,size(W,2));
    % Cluster the subproblems
    G = kmeans(W,K);
    G = arrayfun(@(S)find(G==S),1:K,'UniformOutput',false);

    %% Detect the neighbours of each solution
    B = pdist2(W,W);
    [~,B] = sort(B,2);
    B = B(:,1:T);

    %% Generate random population
    Population = Global.Initialization();
    Z = min(Population.objs,[],1);
    
    %% Initial the CMA model
    sk    = cellfun(@(S)S(randi(length(S))),G);
    xk    = Population(sk).decs;
    Sigma = struct('s',num2cell(sk),'x',num2cell(xk,2)','sigma',0.5,'C',eye(Global.D),'pc',0,'ps',0);

    %% Optimization
    while Global.NotTermination(Population)
        for s = 1 : Global.N
            k = find([Sigma.s]==s);
            if ~isempty(k)
                P = B(s,randperm(size(B,2)));
                % Generate offsprings by CMA
                Offspring = INDIVIDUAL(mvnrnd(Sigma(k).x,Sigma(k).sigma^2*Sigma(k).C,4+floor(3*log(Global.D))));
                % Sort the parent and offsprings
                Combine   = [Offspring,Population(s)];
                [~,rank]  = sort(max(abs(Combine.objs-repmat(Z,length(Combine),1)).*repmat(W(s,:),length(Combine),1),[],2));
                % Update the CMA model
                Sigma(k)  = UpdateCMA(Combine(rank).decs,Sigma(k),Global.gen);
                if isempty(Sigma(k).s)
                    sk = G{k}(randi(length(G{k})));
                    Sigma(k).s = sk;
                    Sigma(k).x = Population(sk).dec;
                end
            else
                % Generate an offspring by DE
                if rand < 0.9
                    P = B(s,randperm(size(B,2)));
                else
                    P = randperm(Global.N);
                end
                Offspring = DE(Population(s),Population(P(1)),Population(P(2)));
            end
            for x = 1 : length(Offspring)
                % Update the ideal point
                Z = min(Z,Offspring(x).obj);
                % Update the solutions in P by Tchebycheff approach
                g_old = max(abs(Population(P).objs-repmat(Z,length(P),1)).*W(P,:),[],2);
                g_new = max(repmat(abs(Offspring(x).obj-Z),length(P),1).*W(P,:),[],2);
                Population(P(find(g_old>=g_new,2))) = Offspring(x);
            end
        end
    end
end