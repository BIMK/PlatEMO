function MOEADDRA(Global)
% <algorithm> <H-N>
% The Performance of a New Version of MOEA/D on CEC09 Unconstrained MOP
% Test Instances
% operator --- DE

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate the weight vectors
    [W,Global.N] = UniformPoint(Global.N,Global.M);
    % Size of neighborhood
    T  = ceil(Global.N/10);
    % Maximum number of solutions replaced by each offspring
    nr = ceil(Global.N/100);

    %% Detect the neighbours of each solution
    B = pdist2(W,W);
    [~,B] = sort(B,2);
    B = B(:,1:T);

    %% Generate random population
    Population = Global.Initialization();
    Z          = min(Population.objs,[],1);
    % Utility for each subproblem
    Pi = ones(Global.N,1);
    % Old Tchebycheff function value of each solution on its subproblem
    oldObj = max(abs((Population.objs-repmat(Z,Global.N,1)).*W),[],2);

    %% Optimization
    while Global.NotTermination(Population)
        for subgeneration = 1 : 5
            % Choose I
            Bounday = find(sum(W<1e-3,2)==Global.M-1)';
            I = [Bounday,TournamentSelection(10,floor(Global.N/5)-length(Bounday),-Pi)];

            % For each solution in I
            for i = I
                % Choose the parents
                if rand < 0.9
                    P = B(i,randperm(size(B,2)));
                else
                    P = randperm(Global.N);
                end

                % Generate an offspring
                Offspring = Global.Variation(Population([i,P(1:2)]),1,@DE);

                % Update the ideal point
                Z = min(Z,Offspring.obj);

                % Update the solutions in P by Tchebycheff approach
                g_old = max(abs(Population(P).objs-repmat(Z,length(P),1)).*W(P,:),[],2);
                g_new = max(repmat(abs(Offspring.obj-Z),length(P),1).*W(P,:),[],2);
                Population(P(find(g_old>=g_new,nr))) = Offspring;
            end
        end
        if ~mod(Global.gen,10)
            % Update Pi for each solution
            newObj    = max(abs((Population.objs-repmat(Z,Global.N,1)).*W),[],2);
            DELTA     = (oldObj-newObj)./oldObj;
            Temp      = DELTA < 0.001;
            Pi(~Temp) = 1;
            Pi(Temp)  = (0.95+0.05*DELTA(Temp)/0.001).*Pi(Temp);
            oldObj    = newObj;
        end
    end
end