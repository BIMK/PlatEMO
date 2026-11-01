classdef CMOEARDA < ALGORITHM
% <2027> <multi> <real> <constrained>
% Constrained multi-objective evolutionary algorithm with role dynamic adjustment
% reward_step ---   5 --- Reward step
% build_step  --- 200 --- DQN model build step
% update_step --- 200 --- DQN model update step

%------------------------------- Reference --------------------------------
% Z. Zeng, Y. Tian, X. Yan, Q. Luo, Q. Dang, and S. Yang. Three-stage
% multi-population constrained multi-objective evolutionary algorithm with
% role dynamic adjustment. Expert Systems with Applications, 2027, 334:
% 134748.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhiqiang Zeng (email: zhiqiang.zeng@outlook.com)

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [reward_step,build_step,update_step] = Algorithm.ParameterSet(5,200,200);
            step        = 0;
            g           = 0;
            num_op      = 2;    % total number of action
            Reb_Data    = [];   % data for DQN
            Sample_Data = [];   % data for DAE
            model_built = 0;    % model build flag
            greedy      = 0.95; % parameter of greedy selection(epsilon value) 
            gama        = 0.95; % discount factor for DNQ
            rate_temp   = [];

            %% Generate random population
            OCPO        = Problem.Initialization(); % objective–constraint balancing population
            UPO         = OCPO;     % unconstrained exploration population
            DEPO        = OCPO;     % diversity enhancement population
            CPO         = OCPO;     % constrained search population
            FPO         = OCPO;     % feasibility exploitation population
            Archive     = OCPO;
            zmin        = min(OCPO.objs,[],1) - 1e-6;
            Sample_Data = [Sample_Data;Archive.objs];
            old_rate    = 0;
            action_op   = randi(num_op);

            %% Optimization
            while Algorithm.NotTerminated(OCPO)
                g    = g + 1;
                step = step + 1;

                if ~model_built && g>build_step         % build model
                    if size(Reb_Data,1) < update_step   % update experience replay buffer
                        use_data = 1:size(Reb_Data,1);
                    else
                        use_data = randperm(size(Reb_Data,1),update_step);
                    end
                    tr_x = Reb_Data(use_data,1:Problem.N+1);
                    tr_y = Reb_Data(use_data,Problem.N+2);
                    operator_net = DQN_Train(tr_x,tr_y);    % build DQN model
                    model_built  = 1;
                end

                if step == 1    % build DAE model
                    dae   = DAEComplete(Problem.M,1);
                    dae.train(Sample_Data);
                    state = dae.reduce(Archive.objs);
                    state = state';
                end

                if Problem.FE/Problem.maxFE < 1/3
                    UPO_selected  = UPO(unidrnd(length(UPO), [1, floor(Problem.N/3)]));
                    DEPO_selected = DEPO(unidrnd(length(DEPO), [1, floor(Problem.N/3)]));

                    MatingPool_Pop1_1 = randperm(length(UPO_selected));
                    MatingPool_Pop1_2 = randperm(length(DEPO_selected));

                    if action_op == 1   % Perform DE search or GA search
                        Offspring1 = DE_Search(UPO_selected,Problem,length(UPO_selected));
                        Offspring2 = DE_Search(DEPO_selected,Problem,length(DEPO_selected));
                        Offspring3 = DE_Search(OCPO,Problem,Problem.N);    
                    else
                        Offspring1 = OperatorGA(Problem, UPO_selected(MatingPool_Pop1_1), {1,20,1,1});
                        Offspring2 = OperatorGA(Problem, DEPO_selected(MatingPool_Pop1_2), {1,20,1,1});
                        Offspring3 = OperatorGA(Problem, OCPO(randperm(length(OCPO))), {1,20,1,1});
                    end

                    Offspring = [Offspring1, Offspring2,Offspring3];
                    zmin      = min(zmin, min(Offspring.objs, [], 1) - 1e-6);
                    OCPO      = OCPO_Selection([OCPO,Offspring],Problem.N,Problem);     % update OCPO
                    [UPO,~]   = UPO_Selection([UPO,Offspring],Problem.N);               % update UPO
                    DEPO      = DEPO_Selection(DEPO, Offspring, zmin, Problem.N);       % update DEPO
                    [Archive,rate_re] = Archive_update(Archive, Offspring,  Problem.N); % update Archive

                    rate_temp = [rate_temp,rate_re];
                    CPO       = UPO;
                else 
                    if Problem.FE/Problem.maxFE < 2/3

                        DEPO_selected     = DEPO(unidrnd(length(DEPO), [1, floor(Problem.N/3)]));
                        CPO_selected      = CPO(unidrnd(length(CPO), [1, floor(Problem.N/3)]));
                        MatingPool_Pop2_1 = randperm(length(DEPO_selected));
                        MatingPool_Pop2_2 = randperm(length(CPO_selected));

                        if action_op == 1       % Perform DE search or GA search
                            Offspring2 = DE_Search(DEPO_selected,Problem,length(DEPO_selected));
                            Offspring3 = DE_Search(CPO_selected,Problem,length(CPO_selected));
                            Offspring4 = DE_Search(OCPO,Problem,Problem.N);
                        else
                            Offspring2 = OperatorGA(Problem, DEPO_selected(MatingPool_Pop2_1), {1,20,1,1});
                            Offspring3 = OperatorGA(Problem, CPO_selected(MatingPool_Pop2_2), {1,20,1,1});
                            Offspring4 = OperatorGA(Problem, OCPO(randperm(length(OCPO))), {1,20,1,1});
                        end

                        Offspring = [Offspring2, Offspring3,Offspring4];
                        zmin      = min(zmin, min(Offspring.objs, [], 1) - 1e-6);
                        DEPO      = DEPO_Selection(DEPO, Offspring, zmin, Problem.N);       % update DEPO
                        [CPO,~]   = CPO_Selection([CPO,Offspring],Problem.N);               % update CPO
                        OCPO      = OCPO_Selection([OCPO,Offspring],Problem.N,Problem);     % update OCPO
                        [Archive,rate_re] = Archive_update(Archive, Offspring,  Problem.N); % update Archive

                        rate_temp = [rate_temp,rate_re];
                        FPO       = CPO;
                    else
                        if action_op == 1   % Perform DE search or GA search
                            Offspring1 = DE_Search([OCPO,FPO,DEPO,Archive],Problem,Problem.N);
                            Offspring2 = DE_Search([OCPO,FPO],Problem,floor(Problem.N/3));
                        else
                            combine_pop1 = [OCPO,FPO,DEPO,Archive];
                            Offspring1   = OperatorGA(Problem, combine_pop1(randperm(length(combine_pop1),Problem.N)), {1,20,1,1});
                            combine_pop2 = [OCPO,FPO];
                            Offspring2   = OperatorGA(Problem, combine_pop2(randperm(length(combine_pop2),floor(Problem.N/3))), {1,20,1,1});
                        end

                        [Archive,rate_re] = Archive_update(Archive, [Offspring1,Offspring2],  Problem.N);   % update Archive
                        OCPO = OCPO_Selection([OCPO,Offspring1,Offspring2,Archive],Problem.N,Problem);      % update OCPO
                        FPO  = FPO_Selection(FPO, [Offspring1,Offspring2], floor(Problem.N/3));             % update FPO

                        combine_pop3 = [OCPO,FPO];
                        zmin = min(zmin, min(combine_pop3.objs, [], 1) - 1e-6);
                        DEPO = DEPO_Selection(DEPO, [OCPO,FPO], zmin, floor(Problem.N/3));  % update DEPO

                        rate_temp = [rate_temp,rate_re];
                    end
                end

                action_op_old = action_op;
                if step >= reward_step
                    state_new = dae.reduce(Archive.objs);
                    state_new = state_new';
                    if model_built  % Using DQN to select a new action
                        reward   = 10*(sum(rate_temp)-old_rate);
                        test_x   = state_new;
                        test_op1 = [test_x,1];
                        test_op1 = test_op1';
                        test_op2 = [test_x,2];
                        test_op2 = test_op2';
                        q_op1    = forward(operator_net,dlarray(test_op1,'CB'));
                        q_op2    = forward(operator_net,dlarray(test_op2,'CB'));

                        if rand > greedy
                            action_op = randi(num_op);
                        else
                            if q_op1 > q_op2
                                action_op = 1;
                            elseif q_op1 < q_op2
                                action_op = 2;
                            else
                                action_op = randi(num_op);
                            end
                        end
                        y= reward + gama*max(double(extractdata(q_op1)),double(extractdata(q_op2)));
                        
                    else
                        y         = 10*(sum(rate_temp)-old_rate);
                        action_op = randi(num_op);
                    end

                    current_record_op = [state action_op_old y state_new];
                    Reb_Data          = [Reb_Data;current_record_op];

                    old_rate  = sum(rate_temp);
                    rate_temp = [];
                    if size(Reb_Data,1) > update_step   % Truncate experience replay pool
                        Reb_Data = Reb_Data(2:size(Reb_Data,1),:);
                    end
                    step = 0;
                end
                Sample_Data = [Sample_Data;Archive.objs];
                if size(Sample_Data,1) > 1000   % Truncate sample data
                    Sample_Data = Sample_Data(101:size(Sample_Data,1),:);
                end

                if model_built && mod(g,update_step) == 0
                    % update operator network
                    if size(Reb_Data,1) <= update_step
                        use_data = 1:size(Reb_Data,1);
                    else
                        use_data = randperm(size(Reb_Data,1),update_step);
                    end
                    tr_x         = Reb_Data(use_data,1:Problem.N+1);
                    tr_y         = Reb_Data(use_data,Problem.N+2);
                    operator_net = DQN_Update(operator_net,tr_x,tr_y);%update DQN model
                end

                if Problem.FE >= Problem.maxFE
                    OCPO = Archive;
                end
            end  
        end
    end
end