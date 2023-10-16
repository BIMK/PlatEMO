function [Population,v] = ECSO(Problem,Population,v,gBest,subswarm_index)
% The particle updating strategy of S-ECSO

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    w     = 0.7968;
    c1    = 1.4962;
    c2    = 1.4962;
    x = Population.decs;
    Obj_fitness = Population.objs;
    %% subswarm
    subswarm1_x = x(subswarm_index(1):subswarm_index(2)-1,:);
    subswarm2_x = x(subswarm_index(2):subswarm_index(3)-1,:);
    subswarm3_x = x(subswarm_index(3):size(x,1),:);

    [Front_rank1,~]= NDSort(Obj_fitness(subswarm_index(1):subswarm_index(2)-1,:),inf);
    subswa = x(find(Front_rank1==1),:);
    lbest(1,:) = subswa(randi(size(subswa,1)),:);
    llBest(subswarm_index(1):subswarm_index(2)-1,:) = repmat(lbest(1,:),size(subswarm1_x,1),1);

    [Front_rank2,~]= NDSort(Obj_fitness(subswarm_index(2):subswarm_index(3)-1,:),inf);
    subswa = x(find(Front_rank2==1)+subswarm_index(2)-1,:);
    lbest(2,:) = subswa(randi(size(subswa,1)),:);
    llBest(subswarm_index(2):subswarm_index(3)-1,:) = repmat(lbest(2,:),size(subswarm2_x,1),1);

    [Front_rank3,~]= NDSort(Obj_fitness(subswarm_index(3):size(x,1),:),inf);
    subswa = x(find(Front_rank3==1)+subswarm_index(3)-1,:);
    lbest(3,:) = subswa(randi(size(subswa,1)),:);
    llBest(subswarm_index(3):size(x,1),:) = repmat(lbest(3,:),size(subswarm3_x,1),1);

    randswarm1 = randperm(size(subswarm1_x,1));
    randswarm2 = randperm(size(subswarm2_x,1)) + size(subswarm1_x,1);
    randswarm3 = randperm(size(subswarm3_x,1)) + size(subswarm1_x,1) + size(subswarm2_x,1);

    sub_pair_size = floor(size(x,1)/3);

    [Front_rank,~] = NDSort(Obj_fitness,inf);
    
    %% Compare
    for i = 1 : sub_pair_size
        randswarm   = [randswarm1(i),randswarm2(i),randswarm3(i)];
        R_randswarm = Front_rank(randswarm);
        [~,rrank]   = min(R_randswarm);
        if rrank == 1
            winner_index = randswarm1(i);
            loser1_index = randswarm2(i);
            loser2_index = randswarm3(i);
        elseif rrank == 2
            winner_index = randswarm2(i);
            loser1_index = randswarm1(i);
            loser2_index = randswarm3(i);
        elseif rrank == 3
            winner_index = randswarm3(i);
            loser1_index = randswarm1(i);
            loser2_index = randswarm2(i);
        end

        %% update
        if rand < 0.5
            v(loser1_index,:) = rand*v(loser1_index,:)+rand*(x(winner_index,:)-x(loser1_index,:));
            x(loser1_index,:) = x(loser1_index,:) + v(loser1_index,:);

            v(loser2_index,:) = w*v(loser2_index,:) + c1*rand*(llBest(loser2_index,:)-x(loser2_index,:))...
                + c2*rand*(gBest-x(loser2_index,:));
            x(loser2_index,:) = x(loser2_index,:) + v(loser2_index,:);
        else
            v(loser2_index,:) = rand*v(loser2_index,:)+rand*(x(winner_index,:)-x(loser2_index,:));
            x(loser2_index,:) = x(loser2_index,:) + v(loser2_index,:);

            v(loser1_index,:) = w*v(loser1_index,:) + c1*rand*(llBest(loser1_index,:)-x(loser1_index,:))...
                + c2*rand*(gBest-x(loser1_index,:));
            x(loser1_index,:) = x(loser1_index,:) + v(loser1_index,:);
        end
    end
    
    %% Restrict the range
    N = size(x,1);
    for irange = 1 : N
        Upper_flag   = Problem.upper<x(irange,:);
        Upper_flag_T = sum(Upper_flag);
        if Upper_flag_T > 0
            Upper_index = find(Upper_flag == 1);
            x(irange,Upper_index) = Problem.upper(Upper_index);
        end

        Low_flag = Problem.lower > x(irange,:);
        Low_flag_T = sum(Low_flag);
        if Low_flag_T > 0
            Low_index = find(Low_flag == 1);
            x(irange,Low_index) = Problem.lower(Low_index);
        end
    end
    
    %% Calculate
    Population = Problem.Evaluation(x);
end