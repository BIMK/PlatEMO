classdef CMOQLMT < ALGORITHM
% <multi> <real> <constrained>
% Constrained multi-objective optimization based on Q-learning and multitasking
% delta --- 0.9 --- The probability of choosing parents locally
% nr    ---   2 --- Maximum number of solutions replaced by each offspring

%------------------------------- Reference --------------------------------
% F. Ming, W. Gong, and L. Gao, Adaptive auxiliary task selection for 
% multitasking-assisted constrained multi-objective optimization [feature],
% IEEE Computational Intelligence Magazine, 2023, 18(2): 18-30.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population{1} = Problem.Initialization();
            Population{2} = Problem.Initialization();
            Population{3} = Problem.Initialization();
            alpha1 = 2./(1+exp(1).^(-Problem.FE*10/Problem.maxFE))-1;
            para = ceil(Problem.maxFE/Problem.N)/2 - ceil(Problem.FE/Problem.N);
            Population{4} = Problem.Initialization();
            Zmin       = min(Population{4}.objs,[],1);
            Fmin       = min(Population{4}(all(Population{4}.cons<=0,2)).objs,[],1);
            Fitness{1} = CalFitness(Population{1}.objs,Population{1}.cons);
            Fitness{3} = CalFitness(Population{3}.objs);
            LengthO{1} = Problem.N/2; LengthO{2} = Problem.N/2; LengthO{3} = Problem.N/2; LengthO{4} = Problem.N/2;
            W          = UniformPoint(Problem.N,Problem.M);
            Ra         = 1;
            
            %% For QL
            num_task = 3;
            Q_Table = zeros(num_task,num_task);
            learning = 0.5;
            gama_ql = 0.9;
            alpha_ql = 0.8;
            greedy_ql = 0.9;
            current_state = randi(num_task);

            %% Optimization
            while Algorithm.NotTerminated(Population{1})
                if Problem.FE < learning * Problem.maxFE
                    % choose an action by greedy strategy based on Q table
                    if rand > greedy_ql || ( Q_Table(current_state,1) == Q_Table(current_state,2) && Q_Table(current_state,1) == Q_Table(current_state,3) )
                        action = randi(num_task);
                    else
                        [~,action] = max(Q_Table(current_state,:));
                    end   
                    
                    % Re-rank
                    Population{1} = Population{1}(randperm(Problem.N));
                    Population{2} = Population{2}(randperm(Problem.N));
                    Lia   = ismember(Population{2}.objs,Population{1}.objs, 'rows');
                    gamma = 1-sum(Lia)/Problem.N;
                    [Population{2},Fitness{2},~] = EnvironmentalSelectionT2(Population{2},Problem.N,alpha1,gamma,para);
                    
                    % Offspring Reproduction
                    MatingPool{1} = TournamentSelection(2,2*LengthO{1},Fitness{1});
                    Offspring{1}  = OperatorGAhalf(Problem,Population{1}(MatingPool{1}));
                    MatingPool{2} = TournamentSelection(2,2*LengthO{2},Fitness{2});
                    Offspring{2}  = OperatorGAhalf(Problem,Population{2}(MatingPool{2}));
                    MatingPool{3} = TournamentSelection(2,2*LengthO{3},Fitness{3});
                    Offspring{3}  = OperatorGAhalf(Problem,Population{3}(MatingPool{3}));
                    Nt     = floor(Ra*Problem.N);
                    if length(Population{1}) > Problem.N-Nt
                        MatingPool{4} = [Population{4}(randsample(Problem.N,Nt)),Population{1}(randsample(Problem.N,Problem.N-Nt))];
                    else
                        MatingPool{4} = Population{4}(randsample(Problem.N,Problem.N));
                    end
                    [Mate1,Mate2,Mate3]        = Neighbor_Pairing_Strategy(MatingPool{4},Zmin);
                    if rand > 0.5
                        Offspring{4}  = OperatorDE(Problem,Mate1,Mate2,Mate3);
                    else
                        Offspring{4}  = OperatorDE(Problem,Mate1,Mate2,Mate3,{0.5,0.5,0.5,0.75});
                    end
                    
                    % determine the transfer rate to update Q table
                    [~,~,Next1] = EnvironmentalSelectionT1([Population{1},Offspring{1},Population{action+1},Offspring{action+1}],Problem.N);
                    succ_rate =  (sum(Next1(length(Population{1})+length(Offspring{1})+1:end))) / (length(Population{action+1})+length(Offspring{action+1}));
                    Q_Table(current_state,action) = Q_Table(current_state,action) + alpha_ql * (succ_rate + gama_ql * (max(Q_Table(action,:))) - Q_Table(current_state,action));
                    
                    current_state = action;
                    
                    % Environmental Selection
                    alpha1 = 2./(1+exp(1).^(-Problem.FE*10/Problem.maxFE)) - 1;
                    para  = ceil(Problem.maxFE/Problem.N)/2 - ceil(Problem.FE/Problem.N);
                    Fmin       = min([Fmin;Offspring{4}(all(Offspring{4}.cons<=0,2)).objs],[],1);
                    Zmin       = min([Zmin;Offspring{4}.objs],[],1);
                    [Population{1},Fitness{1},~] = EnvironmentalSelectionT1([Population{1},Offspring{1},Population{current_state+1},Offspring{current_state+1}],Problem.N);
                    [Population{2},~,~] = EnvironmentalSelectionT2([Population{2},Offspring{2}],Problem.N,alpha1,gamma,para);
                    [Population{3},Fitness{3},~] = EnvironmentalSelectionT3([Population{3},Offspring{3}],Problem.N);
                    [Population{4},~] = ICMA_Update([Population{4},Offspring{4}],Problem.N,W,Zmin,Fmin);
                    
                else
                    % choose an action by greedy strategy based on Q table
                    if rand > greedy_ql || ( Q_Table(current_state,1) == Q_Table(current_state,2) && Q_Table(current_state,1) == Q_Table(current_state,3) )
                        action = randi(num_task);
                    else
                        [~,action] = max(Q_Table(current_state,:));
                    end
                    
                    % Offspring Reproduction
                    MatingPool{1} = TournamentSelection(2,2*LengthO{1},Fitness{1});
                    Offspring{1}  = OperatorGAhalf(Problem,Population{1}(MatingPool{1}));
                    
                    if action == 1
                        % Re-rank
                        Population{1} = Population{1}(randperm(Problem.N));
                        Population{2} = Population{2}(randperm(Problem.N));
                        Lia   = ismember(Population{2}.objs,Population{1}.objs, 'rows');
                        gamma = 1-sum(Lia)/Problem.N;
                        [Population{2},Fitness{2},~] = EnvironmentalSelectionT2(Population{2},Problem.N,alpha1,gamma,para);
                        MatingPool{2} = TournamentSelection(2,2*LengthO{2},Fitness{2});
                        Offspring{2}  = OperatorGAhalf(Problem,Population{2}(MatingPool{2}));
                        
                        % determine the transfer rate to update Q table
                        [~,~,Next1] = EnvironmentalSelectionT1([Population{1},Offspring{1},Population{action+1},Offspring{action+1}],Problem.N);
                        succ_rate =  (sum(Next1(length(Population{1})+length(Offspring{1})+1:end))) / (length(Population{action+1})+length(Offspring{action+1}));
                        Q_Table(current_state,action) = Q_Table(current_state,action) + alpha_ql * (succ_rate + gama_ql * (max(Q_Table(action,:))) - Q_Table(current_state,action));
                        
                        current_state = action;
                        
                        % Environmental Selection
                        alpha1 = 2./(1+exp(1).^(-Problem.FE*10/Problem.maxFE)) - 1;
                        para  = ceil(Problem.maxFE/Problem.N)/2 - ceil(Problem.FE/Problem.N);
                        [Population{1},Fitness{1},~] = EnvironmentalSelectionT1([Population{1},Offspring{1},Population{current_state+1},Offspring{current_state+1}],Problem.N);
                        [Population{2},~,~] = EnvironmentalSelectionT2([Population{2},Offspring{2}],Problem.N,alpha1,gamma,para);
                    elseif action == 2
                        MatingPool{3} = TournamentSelection(2,2*LengthO{3},Fitness{3});
                        Offspring{3}  = OperatorGAhalf(Problem,Population{3}(MatingPool{3}));
                        
                        % determine the transfer rate to update Q table
                        [~,~,Next1] = EnvironmentalSelectionT1([Population{1},Offspring{1},Population{action+1},Offspring{action+1}],Problem.N);
                        succ_rate =  (sum(Next1(length(Population{1})+length(Offspring{1})+1:end))) / (length(Population{action+1})+length(Offspring{action+1}));
                        Q_Table(current_state,action) = Q_Table(current_state,action) + alpha_ql * (succ_rate + gama_ql * (max(Q_Table(action,:))) - Q_Table(current_state,action));
                        
                        current_state = action;
                        
                        % Environmental Selection
                        [Population{1},Fitness{1},~] = EnvironmentalSelectionT1([Population{1},Offspring{1},Population{current_state+1},Offspring{current_state+1}],Problem.N);
                        [Population{3},Fitness{3},~] = EnvironmentalSelectionT3([Population{3},Offspring{3}],Problem.N);
                    else
                        Nt     = floor(Ra*Problem.N);
                        if length(Population{1}) > Problem.N-Nt
                            MatingPool{4} = [Population{4}(randsample(Problem.N,Nt)),Population{1}(randsample(Problem.N,Problem.N-Nt))];
                        else
                            MatingPool{4} = Population{4}(randsample(Problem.N,Problem.N));
                        end
          
                        [Mate1,Mate2,Mate3]        = Neighbor_Pairing_Strategy(MatingPool{4},Zmin);
                        if rand > 0.5
                            Offspring{4}  = OperatorDE(Problem,Mate1,Mate2,Mate3);
                        else
                            Offspring{4}  = OperatorDE(Problem,Mate1,Mate2,Mate3,{0.5,0.5,0.5,0.75});
                        end
                        
                        % determine the transfer rate to update Q table
                        [~,~,Next1] = EnvironmentalSelectionT1([Population{1},Offspring{1},Population{action+1},Offspring{action+1}],Problem.N);
                        succ_rate =  (sum(Next1(length(Population{1})+length(Offspring{1})+1:end))) / (length(Population{action+1})+length(Offspring{action+1}));
                        Q_Table(current_state,action) = Q_Table(current_state,action) + alpha_ql * (succ_rate + gama_ql * (max(Q_Table(action,:))) - Q_Table(current_state,action));
                        
                        current_state = action;
                        
                        % Environmental Selection
                        Fmin       = min([Fmin;Offspring{4}(all(Offspring{4}.cons<=0,2)).objs],[],1);
                        Zmin       = min([Zmin;Offspring{4}.objs],[],1);
                        [Population{1},Fitness{1},~] = EnvironmentalSelectionT1([Population{1},Offspring{1},Population{current_state+1},Offspring{current_state+1}],Problem.N);
                        [Population{4},~] = ICMA_Update([Population{4},Offspring{4}],Problem.N,W,Zmin,Fmin);
                    end
                end
                Ra = 1 - Problem.FE/Problem.maxFE;
            end
        end
    end
end