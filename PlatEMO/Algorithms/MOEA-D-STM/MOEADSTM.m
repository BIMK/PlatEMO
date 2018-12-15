function MOEADSTM(Global)
% <algorithm> <M>
% MOEA/D with stable matching

%------------------------------- Reference --------------------------------
% K. Li, Q. Zhang, S. Kwong, M. Li, and R. Wang, Stable matching-based
% selection in evolutionary multiobjective optimization, IEEE Transactions
% on Evolutionary Computation, 2014, 18(6): 909-923.
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
    % Size of neighborhood
    T  = ceil(Global.N/10);

    %% Detect the neighbours of each solution
    B = pdist2(W,W);
    [~,B] = sort(B,2);
    B = B(:,1:T);

    %% Generate random population
    Population = Global.Initialization();
    z          = min(Population.objs,[],1);
    % Utility for each subproblem
    Pi = ones(Global.N,1);
    % Old Tchebycheff function value of each solution on its subproblem
    oldObj = max(abs((Population.objs-repmat(z,Global.N,1))./W),[],2);

    %% Optimization
    while Global.NotTermination(Population)
        for subgeneration = 1 : 5
            % Choose I
            Bounday = find(sum(W<1e-3,2)==Global.M-1)';
            I = [Bounday,TournamentSelection(10,floor(Global.N/5)-length(Bounday),-Pi)];
            % Generate an offspring for each solution in I
            P = zeros(length(I),3);
            for i = 1 : length(I)
                % Choose the parents
                if rand < 0.9
                    P(i,:) = B(I(i),randperm(size(B,2),3));
                else
                    P(i,:) = randperm(Global.N,3);
                end
            end
            Offspring = DE(Population(P(:,1)),Population(P(:,2)),Population(P(:,3)));
            z         = min([z;Offspring.objs],[],1);
            
            % STM selection
            Population = STM([Population,Offspring],W,z,max(Population.objs,[],1));
        end
        if ~mod(Global.gen,10)
            % Update Pi for each solution
            newObj    = max(abs((Population.objs-repmat(z,Global.N,1))./W),[],2);
            DELTA     = oldObj - newObj;
            Temp      = DELTA < 0.001;
            Pi(~Temp) = 1;
            Pi(Temp)  = (0.95+0.05*DELTA(Temp)/0.001).*Pi(Temp);
            oldObj    = newObj;
        end
    end
end