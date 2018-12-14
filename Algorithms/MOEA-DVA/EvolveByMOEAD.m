function EvolveByMOEAD(Global,Population,W)
% Uniformity optimization by MOEA/D

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Detect the neighbours of each solution
    B = pdist2(W,W);
    [~,B] = sort(B,2);
    B = B(:,1:ceil(Global.N/10));
    
    %% Associate each subproblem with one solution
    % The ideal point
    Z = min(Population.objs,[],1);
    % The value of each solution on each subproblem (modified Tchebycheff approach)
    g = zeros(Global.N);
    for i = 1 : Global.N
        g(i,:) = max(repmat(abs(Population(i).obj-Z),Global.N,1)./W,[],2)';
    end
    [~,rank] = sort(g,2);
    % The index of solution which each subproblem associated with
    associate = zeros(1,Global.N);
    for i = 1 : Global.N
        x = find(~associate(rank(i,:)),1);
        associate(rank(i,x)) = i;
    end
    Population = Population(associate);
    
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
            
            % Update the ideal point
            Z = min(Z,Offspring.obj);

            % Update the solutions in P by modified Tchebycheff approach
            g_old = max(abs(Population(P).objs-repmat(Z,length(P),1))./W(P,:),[],2);
            g_new = max(repmat(abs(Offspring.obj-Z),length(P),1)./W(P,:),[],2);
            Population(P(g_old>=g_new)) = Offspring;
        end
    end
end