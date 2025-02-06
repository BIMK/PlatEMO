classdef DRLOSEMCMO < ALGORITHM
% <2024> <multi> <real/integer/label/binary/permutation> <constrained>
% EMCMO with deep reinforcement learning-assisted operator selection

%------------------------------- Reference --------------------------------
% F. Ming, W. Gong, L. Wang, and Y. Jin. Constrained multi-objective
% optimization with deep reinforcement learning assisted operator
% selection. IEEE/CAA Journal of Automatica Sinica, 2024, 11(4): 919-931.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming (email: 20151000334@cug.edu.cn)

    methods
        function main(Algorithm,Problem)
            Population{1}  = Problem.Initialization();
            Population{2}  = Problem.Initialization();
            Fitness{1}     = CalFitness(Population{1}.objs,Population{1}.cons);
            Fitness{2}     = CalFitness(Population{2}.objs);
            transfer_state = 0;
            cnt            = 0;
            
            %% For DQL
            Data = [];
            num_operator = 2;
            model_built = 0;
            count = 0;
            greedy = 0.95;
            gama = 0.9;
            
            %% Optimization
            while Algorithm.NotTerminated(Population{1})
                gen = ceil(Problem.FE/(2*Problem.N));
                
                % old state
                average_f = sum(sum(Population{1}.objs,2))/length(Population{1});
                cv = overall_cv(Population{1}.cons);
                average_cv = sum(cv)/length(Population{1});
                f_max = max(Population{1}.objs,[],1);
                f_min = min(Population{1}.objs,[],1);
                average_d = sum(f_max-f_min);
                
                if gen <= 200
                    %% Exploring stage, choose operator randomly
                    operator = randi(num_operator);
                else
                    %% Learning stage, choose by Deep Q-net
                    % choose an action based on the trained net
                    if ~model_built
                        % build model here
                        use_data = randperm(length(Data),200);
                        tr_x = Data(use_data,1:4);
                        [tr_xx,ps] = mapminmax(tr_x');tr_xx=tr_xx';
                        tr_y = Data(use_data,5:8);
                        [tr_yy,qs] = mapminmax(tr_y');tr_yy=tr_yy';
                        Params.ps  = ps;Params.qs=qs;
                        [net,Params] = trainmodel(tr_xx,tr_yy,Params);
                        model_built = 1;
                        operator = randi(num_operator);
                    else
                        % use the model to choose action
                        if rand > greedy
                            operator = randi(num_operator);
                        else
                            test_x1 = [average_f,average_cv,average_d,1];
                            test_x2 = [average_f,average_cv,average_d,2];
                            ps=Params.ps;
                            qs=Params.qs;
                            x1=mapminmax('apply',test_x1',ps);x1=x1';
                            x2=mapminmax('apply',test_x2',ps);x2=x2';
                            reward = testNet(x1,net,Params);
                            reward=mapminmax('reverse',reward',qs);reward=reward';
                            succ2 = testNet(x2,net,Params);
                            succ2=mapminmax('reverse',succ2',qs);succ2=succ2';
                            succ = [reward;succ2];
                            [~,operator] = max(succ(:,1));
                        end
                    end
                end
                
                cnt =cnt+1;
                if   transfer_state == 0
                    for i = 1: 2
                        if operator == 1
                            valOffspring{i} = OperatorGAhalf(Problem,Population{i}(randi(Problem.N,1,Problem.N)));
                        else 
                            valOffspring{i} = OperatorDE(Problem,Population{i},Population{i}(randi(ceil(Problem.N),1,Problem.N)),Population{i}(randi(ceil(Problem.N),1,Problem.N)));
                        end
                    end
                    
                    for i = 1:2
                        if i == 1
                            [Population{i},Fitness{i},~] = EnvironmentalSelection( [Population{1},valOffspring{1},valOffspring{2}],Problem.N,i);
                        else
                            [Population{i},Fitness{i},~] = EnvironmentalSelection( [Population{2},valOffspring{2},valOffspring{1}],Problem.N,i);
                        end
                    end
                    
                    if Problem.FE/Problem.maxFE >=0.2
                        transfer_state = 1;
                    end
                    
                else
                    
                    for i = 1: 2
                        if operator == 1
                            MatingPool = TournamentSelection(2,Problem.N,Fitness{i});
                            valOffspring{i} = OperatorGAhalf(Problem,Population{i}(MatingPool));
                        else
                            MatingPool = TournamentSelection(2,2*Problem.N,Fitness{i});
                            valOffspring{i} = OperatorDE(Problem,Population{i},Population{i}(MatingPool(1:end/2)),Population{i}(MatingPool(end/2+1:end)));
                        end
                    end
                    [~,~,Next] = EnvironmentalSelection( [Population{2},valOffspring{2}],Problem.N,1);
                    succ_rate(1,cnt) =  (sum(Next(1:Problem.N))/100) - (sum(Next(Problem.N+1:end))/50);
                    
                    [~,~,Next] = EnvironmentalSelection( [Population{1},valOffspring{1}],Problem.N,2);
                    succ_rate(2,cnt) =  (sum(Next(1:Problem.N))/100) - (sum(Next(Problem.N+1:end))/50);
                    
                    for i = 1:2
                        if   succ_rate(i,cnt) >0
                            rand_number = randperm(Problem.N);
                            [Population{i},Fitness{i},~] = EnvironmentalSelection( [Population{i},valOffspring{i},Population{2/i}(rand_number(1:Problem.N/2))],Problem.N,i);
                        else
                            [Population{i},Fitness{i},~] = EnvironmentalSelection( [Population{i},valOffspring{i},valOffspring{2/i}],Problem.N,i);
                        end
                    end
                end
                
                % new stage
                average_f1 = sum(sum(Population{1}.objs,2))/length(Population{1});
                cv1 = overall_cv(Population{1}.cons);
                average_cv1 = sum(cv1)/length(Population{1});
                f_max1 = max(Population{1}.objs,[],1);
                f_min1 = min(Population{1}.objs,[],1);
                average_d1 = sum(f_max1-f_min1);
                
                %% Update experience replay
                reward = (average_f1 + average_cv1 + average_d1) - (average_f + average_cv + average_d);
                current_record = [average_f average_cv average_d operator reward,average_f1 average_cv1 average_d1];
                Data = [Data;current_record];
                if size(Data,1) > 500
                    Data(end,:) = [];
                end
                
                %% Update Q-net
                % update net model every 50 generations
                if model_built
                    count = count + 1;
                    if count > 50
                        % update model here
                        qs=Params.qs;
                        use_data = randperm(length(Data),200);
                        tr_x = Data(use_data,1:4);
                        [tr_xx,ps] = mapminmax(tr_x');tr_xx=tr_xx';
                        reward = testNet(tr_xx,net,Params);
                        reward=mapminmax('reverse',reward',qs);reward=reward';
                        succ = reward(:,1);
                        tr_yy = Data(use_data,5)+gama*max(succ);
                        [tr_yy,qs] = mapminmax(tr_yy');tr_yy=tr_yy';
                        Params.ps  = ps;Params.qs=qs;
                        net = updatemodel(tr_xx,tr_yy,Params,net);
                        count = 0;
                    end
                end
            end
        end
    end
end

function result = overall_cv(cv)
% The overall constraint violation

    cv(cv <= 0) = 0;
    cv = abs(cv);
    result = sum(cv,2);
end