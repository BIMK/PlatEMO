function MOEADAWA(Global)
% <algorithm> <M>
% MOEA/D with adaptive weight adjustment
% rate_update_weight --- 0.05 --- Ratio of updated weight vectors
% rate_evol          ---  0.8 --- Ratio of iterations to evolve with only MOEA/D
% wag                ---  100 --- Iteration interval of utilizing AWA

%------------------------------- Reference --------------------------------
% Y. Qi, X. Ma, F. Liu, L. Jiao, J. Sun, and J. Wu, MOEA/D with adaptive
% weight adjustment, Evolutionary Computation, 2014, 22(2): 231-264.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [rate_update_weight,rate_evol,wag] = Global.ParameterSet(0.05,0.8,100);

    %% Generate the weight vectors
    [W,Global.N] = UniformPoint(Global.N,Global.M);
    % Transformation on W
    W = 1./W./repmat(sum(1./W,2),1,size(W,2));
    % Size of neighborhood
    T = ceil(Global.N/10);
    % Maximum number of solutions replaced by each offspring
    nr = ceil(Global.N/100);
    % Size of external elite
    nEP = ceil(Global.N*1.5);

    %% Detect the neighbours of each solution
    B = pdist2(W,W);
    [~,B] = sort(B,2);
    B = B(:,1:T);

    %% Generate random population
    Population = Global.Initialization();
    Z          = min(Population.objs,[],1);
    Pi         = ones(Global.N,1);
    oldObj     = max(abs((Population.objs-repmat(Z,Global.N,1)).*W),[],2);

    %% Optimization
    EP = [];
    while Global.NotTermination(Population)
        if ~mod(Global.gen,10)
            % Allocation of computing resources
            newObj    = max(abs((Population.objs-repmat(Z,Global.N,1)).*W),[],2);
            DELTA     = (oldObj-newObj)./oldObj;
            Temp      = DELTA <= 0.001;
            Pi(~Temp) = 1;
            Pi(Temp)  = (0.95+0.05*DELTA(Temp)/0.001).*Pi(Temp);
            oldObj    = newObj;
        end
        for subgeneration = 1 : 5
            % Choose I
            Bounday = find(sum(W<1e-3,2)==1)';
            I = [Bounday,TournamentSelection(10,floor(Global.N/5)-length(Bounday),-Pi)];

            % Evolve each solution in I
            Offsprings(1:length(I)) = INDIVIDUAL();
            for i = 1 : length(I)
                % Choose the parents
                if rand < 0.9
                    P = B(I(i),randperm(size(B,2)));
                else
                    P = randperm(Global.N);
                end

                % Generate an offspring
                Offsprings(i) = GAhalf(Population(P(1:2)));

                % Update the ideal point
                Z = min(Z,Offsprings(i).obj);

                % Update the solutions in P by Tchebycheff approach
                g_old = max(abs(Population(P).objs-repmat(Z,length(P),1)).*W(P,:),[],2);
                g_new = max(repmat(abs(Offsprings(i).obj-Z),length(P),1).*W(P,:),[],2);
                Population(P(find(g_old>=g_new,nr))) = Offsprings(i);
            end
        end
        if Global.gen >= rate_evol*Global.maxgen
            % Adaptive weight adjustment
            if isempty(EP)
                EP = updateEP(Population,Offsprings,nEP);
            else
                EP = updateEP(EP,Offsprings,nEP);
            end
            if ~mod(Global.gen,wag/5)
                [Population,W] = updateWeight(Population,W,Z,EP,rate_update_weight*Global.N);
            end
        end
    end
end