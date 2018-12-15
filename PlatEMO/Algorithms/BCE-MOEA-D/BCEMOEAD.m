function BCEMOEAD(Global)
% <algorithm> <B>
% Bi-criterion evolution based MOEA/D

%------------------------------- Reference --------------------------------
% M. Li, S. Yang, and X. Liu, Pareto or non-Pareto: Bi-criterion evolution
% in multiobjective optimization, IEEE Transactions on Evolutionary
% Computation, 2016, 20(5): 645-665.
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
    T  = ceil(Global.N/10);
    nr = ceil(Global.N/100);

    %% Detect the neighbours of each solution
    B = pdist2(W,W);
    [~,B] = sort(B,2);
    B = B(:,1:T);
    
    %% Generate random population
    NPC = Global.Initialization(size(W,1));
    Z   = min(NPC.objs,[],1);
    [PC,nND] = PCSelection(NPC,Global.N);

    %% Optimization
    while Global.NotTermination(PC)
        % PC evolving
        NewPC = Exploration(PC,NPC,nND,Global.N);
        
        % NPC selection
        for i = 1 : length(NewPC)
            % Update the ideal point
            Z = min(Z,NewPC(i).obj);
            % Update at most one solution in NPC
            P     = randperm(length(NPC));
            g_old = max(abs(NPC(P).objs-repmat(Z,length(P),1))./W(P,:),[],2);
            g_new = max(repmat(abs(NewPC(i).obj-Z),length(P),1)./W(P,:),[],2);
            NPC(P(find(g_old>=g_new,1))) = NewPC(i);
        end
        
        % NPC evolving
        NewNPC(1:length(NPC)) = INDIVIDUAL();
        for i = 1 : length(NPC)
            % Choose the parents
            if rand < 0.9
                P = B(i,randperm(size(B,2)));
            else
                P = randperm(length(NPC));
            end
            % Generate an offspring
            NewNPC(i) = GAhalf(NPC(P(1:2)));
            % Update the ideal point
            Z = min(Z,NewNPC(i).obj);
            % Update the solutions in P by modified Tchebycheff approach
            g_old = max(abs(NPC(P).objs-repmat(Z,length(P),1))./W(P,:),[],2);
            g_new = max(repmat(abs(NewNPC(i).obj-Z),length(P),1)./W(P,:),[],2);
            NPC(P(find(g_old>=g_new,nr))) = NewNPC(i);
        end
     
        % PC selection
        [PC,nND] = PCSelection([PC,NewNPC,NewPC],Global.N);
    end
end