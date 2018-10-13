function MOEADDU(Global)
% <algorithm> <H-N>
% Balancing Convergence and Diversity in Decomposition-Based Many-Objective
% Optimizers
% delta --- 0.9 --- The probability of choosing parents locally
% K     ---   5 --- Number of nearest weight vectors

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [delta,K] = Global.ParameterSet(0.9,5);

    %% Generate the weight vectors
    [W,Global.N] = UniformPoint(Global.N,Global.M);
    T = ceil(Global.N/10);

    %% Detect the neighbours of each solution
    B = pdist2(W,W);
    [~,B] = sort(B,2);
    B = B(:,1:T);

    %% Generate random population
    Population = Global.Initialization();
    % Ideal and nadir points
    z    = min(Population.objs,[],1);
    znad = max(Population.objs,[],1);

    %% Optimization
    while Global.NotTermination(Population)
        % Normalization
        [~,z,znad] = Normalization(Population.objs,z,znad);
        % For each solution
        for i = 1 : Global.N
            % Choose the parent
            if rand < delta
                P = B(i,randi(size(B,2)));
            else
                P = randi(Global.N);
            end

            % Generate an offspring
            Offspring = Global.Variation(Population([i,P]),1);

            % Find the K nearest weight vectors of the offspring
            [~,rank] = sort(1-pdist2(Offspring.obj,W,'cosine'),'descend');
            P        = rank(1:K);

            % Update the K nearest parents by modified Tchebycheff approach
            g_old = max(abs(Population(P).objs-repmat(z,length(P),1))./repmat(znad-z,K,1)./W(P,:),[],2);
            g_new = max(repmat(abs(Offspring.obj-z)./(znad-z),length(P),1)./W(P,:),[],2);
            Population(P(find(g_old>=g_new,1))) = Offspring;
        end
    end
end