function [SwarmRBF,DeltaRBF] = RBFOperator(net,Demons,SwarmRBF,DeltaRBF,Gbest,Problem)
% OperatorPSO - The operator of particle swarm optimization.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    factor = 0;
    [N,D]  = size(SwarmRBF);
    D      = D - 1;
    BU     = Problem.upper;
    BD     = Problem.lower;
    
    %% SL-PSO
    for i = 1 : N
        % Choose an individual to learn one dimension to follow
        BetterInd = find(Demons(:,D+1)<SwarmRBF(i,D+1));
        ChoseId   = zeros(1,D);
        if length(BetterInd) == 1
            ChoseId = Demons(BetterInd(1),1:D);
        elseif length(BetterInd) == 2
            Mix = rand(1,D) > 0.5;
            ChoseId(Mix) = Demons(BetterInd(1),Mix);
            ChoseId(~Mix) = Demons(BetterInd(2),~Mix);
        elseif length(BetterInd) > 2
            for t = 1 : D
                ChoseId(t) = Demons(BetterInd(randperm(length(BetterInd),1)),t);
            end
        end
        if ~isempty(BetterInd)
            DeltaRBF(i,:) = rand(1,D).*DeltaRBF(i,:) + rand(1,D).*(ChoseId-SwarmRBF(i,1:D));
            % Update
            DeltaRBF(i,1:D) = max(DeltaRBF(i,1:D),BD);
            DeltaRBF(i,1:D) = min(DeltaRBF(i,1:D),BU);
            SwarmRBF(i,1:D) = DeltaRBF(i,1:D);
        end
    end
    %% Repair
    SwarmRBF(:,1:D) = max(SwarmRBF(:,1:D),repmat(Problem.lower,N,1));
    SwarmRBF(:,1:D) = min(SwarmRBF(:,1:D),repmat(Problem.upper,N,1));
    % Approximate the fitness of each particle
    SwarmRBF(:,1+D) = sim(net,SwarmRBF(:,1:D)')';
end