function MOEADPaS(Global)
% <algorithm> <M>
% MOEA/D with Pareto adaptive scalarizing approximation

%------------------------------- Reference --------------------------------
% R. Wang, Q. Zhang, and T. Zhang, Decomposition-based algorithms using
% Pareto adaptive scalarizing methods, IEEE Transactions on Evolutionary
% Computation, 2016, 20(6): 821-837.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate the weight vectors
    [W,Global.N] = UniformPoint(Global.N,Global.M);
    
    %% Initialize the norm of each weight vector
    p = ones(Global.N,1);

    %% Detect the neighbours of each solution
    T = ceil(Global.N/10);
    B = pdist2(W,W);
    [~,B] = sort(B,2);
    B = B(:,1:T);

    %% Generate random population
    Population = Global.Initialization();
    % Ideal and nadir points
    z    = min(Population.objs,[],1);
    znad = max(Population(NDSort(Population.objs,1)==1).objs,[],1);

    %% Optimization
    while Global.NotTermination(Population)
        % For each solution
        for i = 1 : Global.N
            % Choose the parents
            if rand < 0.9
                P = B(i,randperm(size(B,2)));
            else
                P = randperm(Global.N);
            end

            % Generate an offspring
            Offspring = DE(Population(i),Population(P(1)),Population(P(2)));

            % Update the solutions in P by the scalarizing method
            INF   = p(P) == inf;
            g_old = zeros(1,length(P));
            g_new = zeros(1,length(P));
            if any(INF)
                g_old(INF) = max((Population(P(INF)).objs-repmat(z,sum(INF),1))./repmat(znad-z,sum(INF),1)./W(P(INF),:),[],2);
                g_new(INF) = max(repmat((Offspring.obj-z)./(znad-z),sum(INF),1)./W(P(INF),:),[],2);
            end
            if any(~INF)
                g_old(~INF) = sum(((Population(P(~INF)).objs-repmat(z,sum(~INF),1))./repmat(znad-z,sum(~INF),1)./W(P(~INF),:)).^repmat(p(P(~INF)),1,Global.M),2).^(1./p(P(~INF)));
                g_new(~INF) = sum((repmat((Offspring.obj-z)./(znad-z),sum(~INF),1)./W(P(~INF),:)).^repmat(p(P(~INF)),1,Global.M),2).^(1./p(P(~INF)));
            end
            Population(P(find(g_old>g_new,ceil(0.1*T)))) = Offspring;
        end
        
        % Update the ideal point and nadir point
        z    = min([z;Population.objs],[],1);
        znad = max(Population(NDSort(Population.objs,1)==1).objs,[],1);
        
        % Update the norm of each weight vector
        P    = [1,2,3,4,5,6,7,8,9,10,inf];
        nObj = (Population.objs-repmat(z,Global.N,1))./repmat(znad-z,Global.N,1);
        for i = find(rand(1,Global.N)>=Global.gen/Global.maxgen)
            g = zeros(Global.N,length(P));
            for j = 1 : length(P)
                if P(j) ~= inf
                    g(:,j) = sum((nObj./repmat(W(i,:),Global.N,1)).^P(j),2).^(1./P(j));
                else
                    g(:,j) = max((nObj./repmat(W(i,:),Global.N,1)),[],2);
                end
            end
            [~,ZK] = min(g,[],1);
            [~,Z]  = min(sqrt(1-(1-pdist2(nObj(ZK,:),W(i,:),'cosine')).^2).*sqrt(sum(nObj(ZK,:).^2,2)));
            p(i)   = P(Z);
        end
    end
end