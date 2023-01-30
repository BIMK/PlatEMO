function [Population,Position,Velocity,Gbest,Pbest] = UpdatePosition(Problem,Population,Position,Velocity,Pbest,Gbest,currFES,maxFES,Wnc,Pr,model,eta)
% Update Position for L and S swarm in SAMSO

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Determine the size of each swarm
    [N,D] = size(Velocity);
    N2    = min(N,2+ceil(N*((maxFES-currFES)/maxFES)^Wnc));
    N1    = N - N2;
     
    %% Produce offspring for L and S swarm
    tVelocity = Velocity;
    tPosition = Position;
    % S swarm
    if N1 ~= 0
        IDX = 1:N1;

        % Standard PSO
        label = rand(N1,1) < Pr;
        w  = 0.792;
        c1 = 1.491;
        c2 = 1.491;
        r1 = repmat(rand(sum(label),1),1,D);
        r2 = repmat(rand(sum(label),1),1,D);
        tVelocity(IDX(label),:) = w*Velocity(IDX(label),:) + c1*r1.*(Pbest(IDX(label),1:D)-Position(IDX(label),1:D))+...
            c2*r2.*(repmat(Gbest(1:D),sum(label),1)-Position(IDX(label),1:D));
        tPosition(IDX(label),1:D) = Position(IDX(label),1:D) + tVelocity(IDX(label),:);
        
        % Eigencoordinate System
        [~,I]  = sort(Population.objs,'ascend');
        num    = min(length(Population.objs),2*D);
        TopDec = Population(I(1:num)).decs;
        B      = orth(cov(TopDec));
        remain = find(~label);
        if ~isempty(remain)
            num = length(remain);
            R1  = rand(1,num);
            R2  = rand(1,num);
            for i = 1 : num
                    tVelocity(IDX(remain(i)),:) = w*Velocity(IDX(remain(i)),:) + c1*(Pbest(IDX(remain(i)),1:D)-Position(IDX(remain(i)),1:D))*B*R1(i)*B'+...
                c1*(Gbest(1:D)-Position(IDX(remain(i)),1:D))*B*R2(i)*B';
            end
        end
    end
    % L-Swarm
    learner1 = 1:1:N2;
    learner2 = randperm(N2);
    while any(learner1==learner2)
        replace  = learner1==learner2;
        learner2(replace) = randperm(N2,sum(replace));
    end
    delta  = Position(N1+learner1,1:D) - Position(N1+learner2,1:D);
    change = Position(N1+learner1,D+1)>=Position(N1+learner2,D+1);
    delta(change) = -delta(change);
    tPosition(N1+1:end,1:D) = Position(N1+1:end,1:D) + rand(N2,1).*delta;
    
    %% Select particles
    % Predict the fitness value
    srgtObj    = rbf_predict(model,Population(1:model.n).decs,tPosition(:,1:D));
    dist       = min(pdist2(tPosition(:,1:D),Population.decs),[],2);
    evaluation = find(dist > eta & srgtObj < Position(:,D+1));
    if ~isempty(evaluation)
        news   = Problem.Evaluation(tPosition(evaluation,1:D));
        tPosition(evaluation,D+1) = news.objs;
        % Update particle
        update = tPosition(evaluation,D+1) < Position(evaluation,D+1);
        Position(evaluation(update),:) = tPosition(evaluation(update),:);
        Velocity(evaluation(update),:) = tVelocity(evaluation(update),:);
        Pbest(evaluation(update),:)    = tPosition(evaluation(update),:);
        if Pbest(evaluation(update),D+1) < Gbest(D+1)
            tPbest   = Pbest(evaluation(update),:);
            [~,best] = min(tPbest(D+1));
            Gbest    = tPbest(best,:);
        end
        Population = [Population,news];
    else
        % No new position is exactly evaluated
        [~,best] = min(tPosition(:,D+1));
        new      = Problem.Evaluation(tPosition(best,1:D));
        Population = [Population,new];
        if new.objs < Position(best,D+1)
            Position(best,1:D) = tPosition(best,1:D);
            Position(best,D+1) = new.objs;
            Velocity(best,:)   = tVelocity(best,:);
            Pbest(best,:)      = Position(best,:);
            if new.objs < Gbest(D+1)
                Gbest = Position(best,:);
            end
        end
    end
end