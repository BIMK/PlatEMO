classdef CMODRL < ALGORITHM
% <2025> <multi> <real/integer/label/binary/permutation> <constrained>
% Constrained multiobjective optimization via deep reinforcement learning
% reward_step --- 5 --- Reward step setting

%------------------------------- Reference --------------------------------
% F. Ming, W. Gong, B. Xue, M. Zhang, and Y. Jin. Automated configuration
% of evolutionary algorithms via deep reinforcement learning for
% constrained multiobjective optimization. IEEE Transactions on
% Cybernetics, 2025, 55(12): 5877-5890.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            reward_step = Algorithm.ParameterSet(5);
            g           = 0;
            Data        = [];
            op_Data     = [];
            PopData     = [];
            model_built = 0;
            greedy      = 0.95;
            gama        = 0.95;
            step        = 0;
            accumulate  = 0;
            num_op      = 2;
            
            %% Generate random population
            Population1 = Problem.Initialization();
            Population2 = Problem.Initialization(ceil(Problem.N/2));
            Archive     = EnvironmentalSelection([Population1,Population2],Problem.N,true);
            initial_cv  = sum(sum(max(0,Population1.cons),2))/Problem.N;
            F_e         = rand/5 * initial_cv;
            
            %% Evaluate solutions based on static environment at initial iteration
            Fitness1 = CalFitness(Population1.objs,Population1.cons);
            Fitness2 = CalFitness(Population2.objs);

            %% Optimization
            while Algorithm.NotTerminated(Population1)
                g = g + 1;
                if g <= 100
                    % Generate offspring
                    MatingPool1 = TournamentSelection(2,Problem.N,Fitness1);
                    MatingPool2 = TournamentSelection(2,Problem.N/2,Fitness2);
                    Offspring1  = OperatorGA(Problem,Population1(MatingPool1));
                    Offspring2  = OperatorGA(Problem,Population2(MatingPool2));
                    
                    % Update archive
                    [Archive,~,~] = EnvironmentalSelection([Archive,Offspring1,Offspring2],Problem.N,true);
                    
                    % Select based on learned environment and static environment
                    [Population1,Fitness1]   = EnvironmentalSelectionFe([Population1,Offspring1,Offspring2],Problem.N,F_e);
                    [Population2,Fitness2,~] = EnvironmentalSelection([Population2,Offspring1,Offspring2],Problem.N/2,false); 
                    
                % Randomly sample in the early 600 generations for exploration
                elseif g <= 600
                    step = step + 1;
                    if step == 1
                        accumulate = 0;
                        
                        % Random F_e and operator
                        average_CV = sum(sum(max(0,Population1.cons),2))/Problem.N;
                        F_e = rand/5 * average_CV;
                        if F_e <= 0
                            F_e = rand/5 * initial_cv * (1-g/ceil(Problem.maxFE/2*Problem.N));
                        end
                        action_op = randi(num_op);
                        
                        % old state
                        % use archive solutions and population1 solutions
                        reward_f  = sum(sum(Archive.objs,2))/Problem.N;
                        reward_cv = sum(sum(max(0,Archive.cons),2))/length(Archive);
                        f_max     = max(Archive.objs,[],1);
                        f_min     = min(Archive.objs,[],1);
                        reward_d  = sum(f_max-f_min);
                        
                        % use population1 solutions
                        % build the DAE/RBM model to estimate the state
                        dae   = DAE(Problem.M,1,10,Problem.N,0.5,0.5,0.1);
                        dae.train(PopData);
                        state = dae.reduce(Population1.objs);
                        state = state';
                        
                    end
                    
                    % Generate offspring
                    if action_op == 1
                        % use ordinary DE
                        MatingPool1 = TournamentSelection(2,2*Problem.N,Fitness1);
                        MatingPool2 = TournamentSelection(2,Problem.N,Fitness2);
                        Offspring1  = OperatorDE(Problem,Population1,Population1(MatingPool1(1:end/2)),Population1(MatingPool1(end/2+1:end)));
                        Offspring2  = OperatorDE(Problem,Population2,Population2(MatingPool2(1:end/2)),Population2(MatingPool2(end/2+1:end)));
                    else
                        % use GA
                        MatingPool1 = TournamentSelection(2,Problem.N,Fitness1);
                        MatingPool2 = TournamentSelection(2,Problem.N/2,Fitness2);
                        Offspring1  = OperatorGA(Problem,Population1(MatingPool1));
                        Offspring2  = OperatorGA(Problem,Population2(MatingPool2));
                    end
                    
                    % Update archive
                    [Archive,~,Next] = EnvironmentalSelection([Archive,Offspring1,Offspring2],Problem.N,true);
                    accumulate = accumulate + sum(Next(length(Archive)+1:length(Archive)+length(Offspring1)))/length(Archive);
                    
                    % Select based on learned environment and static environment
                    [Population1,Fitness1]   = EnvironmentalSelectionFe([Population1,Offspring1,Offspring2],Problem.N,F_e);
                    [Population2,Fitness2,~] = EnvironmentalSelection([Population2,Offspring1,Offspring2],Problem.N/2,false);
                    
                    if step >= reward_step
                        % new state
                        % use population1 solutions
                        state_new = dae.reduce(Population1.objs);
                        state_new = state_new';
                        % use improvement on archive
                        reward_f1  = sum(sum(Archive.objs,2))/Problem.N;
                        reward_cv1 = sum(sum(max(0,Archive.cons),2))/length(Archive);
                        f_max      = max(Archive.objs,[],1);
                        f_min      = min(Archive.objs,[],1);
                        reward_d1  = sum(f_max-f_min);
                        
                        y = 1000 * ((reward_f + reward_cv + reward_d1) - (reward_f1 + reward_cv1 + reward_d) + accumulate/reward_step);
                        
                        % Update records and replay buffer
                        current_record    = [state F_e y state_new];
                        Data              = [Data;current_record];
                        current_record_op = [state action_op y state_new];
                        op_Data           = [op_Data;current_record_op];
                        step              = 0;
                    end
                % use DRL model every 5 generations, means keep F_e and action_op unchanged in 5 generations
                else
                    step = step + 1;
                    if step == 1
                        accumulate = 0;
                        
                        % reward
                        reward_f  = sum(sum(Archive.objs,2))/Problem.N;
                        reward_cv = sum(sum(max(0,Archive.cons),2))/length(Archive);
                        f_max     = max(Archive.objs,[],1);
                        f_min     = min(Archive.objs,[],1);
                        reward_d  = sum(f_max-f_min);
                        
                        % state
                        % build the DAE/RBM model to estimate the state
                        dae   = DAE(Problem.M,1,10,Problem.N,0.5,0.5,0.1);
                        dae.train(PopData);
                        state = dae.reduce(Population1.objs);
                        state = state';
                    end
                    if ~model_built
                        if size(Data,1) < 120
                            use_data = 1:size(Data,1);
                        else
                            use_data = randperm(size(Data,1),120);
                        end
                        % build initial Q-net, also named critic net
                        tr_x       = Data(use_data,1:Problem.N+1);
                        tr_y       = Data(use_data,Problem.N+2);
                        critic_net = TrainCritic(tr_x,tr_y);
                        tr_x       = Data(use_data,1:Problem.N);
                        actor_net  = TrainActor(critic_net,tr_x);
                        
                        % build initial operator network
                        if size(Data,1) < 120
                            use_data = 1:size(op_Data,1);
                        else
                            use_data = randperm(size(op_Data,1),120);
                        end
                        tr_x = op_Data(use_data,1:Problem.N+1);
                        tr_y = op_Data(use_data,Problem.N+2);
                        operator_net = TrainDQN(tr_x,tr_y);
                        
                        % determine actions based on DRL
                        test_x = state;
                        % determine F_e
                        action_es = forward(actor_net,dlarray(test_x','CB'));
                        F_e       = double(extractdata(action_es));
                        
                        % determine operator
                        test_x   = state;
                        test_op1 = [test_x,1];
                        test_op1 = test_op1';
                        test_op2 = [test_x,2];
                        test_op2 = test_op2';
                        q_op1    = forward(operator_net,dlarray(test_op1,'CB'));
                        q_op2    = forward(operator_net,dlarray(test_op2,'CB'));
                        if q_op1 > q_op2
                            action_op = 1;
                        elseif q_op1 < q_op2
                            action_op = 2;
                        else
                            action_op = randi(num_op);
                        end
                        model_built = 1;
                        
                        % Generate offspring
                        if action_op == 1
                            % use ordinary DE
                            MatingPool1 = TournamentSelection(2,2*Problem.N,Fitness1);
                            MatingPool2 = TournamentSelection(2,Problem.N,Fitness2);
                            Offspring1  = OperatorDE(Problem,Population1,Population1(MatingPool1(1:end/2)),Population1(MatingPool1(end/2+1:end)));
                            Offspring2  = OperatorDE(Problem,Population2,Population2(MatingPool2(1:end/2)),Population2(MatingPool2(end/2+1:end)));
                        else
                            % use GA
                            MatingPool1 = TournamentSelection(2,Problem.N,Fitness1);
                            MatingPool2 = TournamentSelection(2,Problem.N/2,Fitness2);
                            Offspring1  = OperatorGA(Problem,Population1(MatingPool1));
                            Offspring2  = OperatorGA(Problem,Population2(MatingPool2));
                        end
                        
                        % Update archive
                        [Archive,~,Next] = EnvironmentalSelection([Archive,Offspring1,Offspring2],Problem.N,true);
                        accumulate       = accumulate + sum(Next(length(Archive)+1:length(Archive)+length(Offspring1)))/length(Archive);
                        
                        % Select based on learned environment and static environment
                        [Population1,Fitness1]   = EnvironmentalSelectionFe([Population1,Offspring1,Offspring2],Problem.N,F_e);
                        [Population2,Fitness2,~] = EnvironmentalSelection([Population2,Offspring1,Offspring2],Problem.N/2,false);
                    else
                        % Generate offspring
                        if action_op == 1
                            % use ordinary DE
                            MatingPool1 = TournamentSelection(2,2*Problem.N,Fitness1);
                            MatingPool2 = TournamentSelection(2,Problem.N,Fitness2);
                            Offspring1  = OperatorDE(Problem,Population1,Population1(MatingPool1(1:end/2)),Population1(MatingPool1(end/2+1:end)));
                            Offspring2  = OperatorDE(Problem,Population2,Population2(MatingPool2(1:end/2)),Population2(MatingPool2(end/2+1:end)));
                        else
                            % use GA
                            MatingPool1 = TournamentSelection(2,Problem.N,Fitness1);
                            MatingPool2 = TournamentSelection(2,Problem.N/2,Fitness2);
                            Offspring1  = OperatorGA(Problem,Population1(MatingPool1));
                            Offspring2  = OperatorGA(Problem,Population2(MatingPool2));
                        end
                        
                        % Update archive
                        [Archive,~,Next] = EnvironmentalSelection([Archive,Offspring1,Offspring2],Problem.N,true);
                        accumulate       = accumulate + sum(Next(length(Archive)+1:length(Archive)+length(Offspring1)))/length(Archive);
                        
                        % Select based on learned environment and static environment
                        [Population1,Fitness1]   = EnvironmentalSelectionFe([Population1,Offspring1,Offspring2],Problem.N,F_e);
                        [Population2,Fitness2,~] = EnvironmentalSelection([Population2,Offspring1,Offspring2],Problem.N/2,false);
                        
                        if step >= reward_step
                            % Update replay buffer
                            % get new state
                            state_new = dae.reduce(Population1.objs);
                            state_new = state_new';
                            
                            % real reward
                            reward_f1  = sum(sum(Archive.objs,2))/Problem.N;
                            reward_cv1 = sum(sum(max(0,Archive.cons),2))/length(Archive);
                            f_max      = max(Archive.objs,[],1);
                            f_min      = min(Archive.objs,[],1);
                            reward_d1  = sum(f_max-f_min);
                            reward     = 1000 * ((reward_f + reward_cv + reward_d1) - (reward_f1 + reward_cv1 + reward_d) + accumulate/reward_step);
                            
                            % estimated reward of actor-critic
                            test_x_actor   = state_new;
                            target_action  = forward(actor_net,dlarray(test_x_actor','CB'));
                            test_x_critic  = [state_new,double(extractdata(target_action))];
                            target_reward  = forward(critic_net,dlarray(test_x_critic','CB'));
                            target_reward  = double(extractdata(target_reward));
                            y              = reward + gama*target_reward;
                            current_record = [state F_e y state_new];
                            Data           = [Data;current_record];
                            if size(Data,1) > 200
                                Data = Data(2:size(Data,1),:);
                            end
                            
                            % estimated reward of operator-net
                            test_x_op1        = [state_new,action_op];
                            target_reward     = forward(operator_net,dlarray(test_x_op1','CB'));
                            target_reward     = double(extractdata(target_reward));
                            y                 = reward + gama*target_reward;
                            current_record_op = [state action_op y state_new];
                            op_Data           = [op_Data;current_record_op];
                            if size(op_Data,1) > 200
                                op_Data = op_Data(2:size(op_Data,1),:);
                            end
                            
                            % Determine new actions
                            if rand > greedy
                                average_CV = sum(sum(max(0,Population1.cons),2))/Problem.N;
                                F_e        = rand/5 * average_CV;
                                if F_e <= 0
                                    F_e = rand/5 * initial_cv * (1-g/ceil(Problem.maxFE/2*Problem.N));
                                end
                                action_op = randi(num_op);
                            else
                                % determine F_e
                                test_x    = state_new;
                                action_es = forward(actor_net,dlarray(test_x','CB'));
                                F_e       = double(extractdata(action_es));
                                
                                % determine operator
                                test_x   = state_new;
                                test_op1 = [test_x,1];
                                test_op1 = test_op1';
                                test_op2 = [test_x,2];
                                test_op2 = test_op2';
                                q_op1    = forward(operator_net,dlarray(test_op1,'CB'));
                                q_op2    = forward(operator_net,dlarray(test_op2,'CB'));
                                if q_op1 > q_op2
                                    action_op = 1;
                                elseif q_op1 < q_op2
                                    action_op = 2;
                                else
                                    action_op = randi(num_op);
                                end
                            end
                            step = 0;
                        end

                    end
                end
                PopData = [PopData;Population1.objs];
                if size(PopData,1) >= 1000
                    PopData = PopData(101:size(PopData,1),:);
                end
                
                % Update networks every 200 generations
                if model_built && mod(g,200) == 0
                    if size(Data,1) < 120
                        use_data = 1:size(Data,1);
                    else
                        use_data = randperm(size(Data,1),120);
                    end
                    % update Q-net
                    tr_x       = Data(use_data,1:Problem.N+1);
                    tr_y       = Data(use_data,Problem.N+2);
                    critic_net = UpdateCritic(critic_net,tr_x,tr_y);
                    % update policy net
                    tr_x      = Data(use_data,1:Problem.N);
                    actor_net = UpdateActor(actor_net,critic_net,tr_x);
                    
                    % update operator network
                    if size(op_Data,1) <= 120
                        use_data = 1:size(op_Data,1);
                    else
                        use_data = randperm(size(op_Data,1),120);
                    end
                    tr_x         = op_Data(use_data,1:Problem.N+1);
                    tr_y         = op_Data(use_data,Problem.N+2);
                    operator_net = UpdateDQN(operator_net,tr_x,tr_y);
                end
                if Problem.FE >= Problem.maxFE
                    Population1 = Archive;
                end
            end
        end
    end
end