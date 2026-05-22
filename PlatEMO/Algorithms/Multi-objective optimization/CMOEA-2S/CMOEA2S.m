classdef CMOEA2S < ALGORITHM
% <2025> <multi> <real> <constrained>
% Constrained MOEA with two types of evolution stages

%------------------------------- Reference --------------------------------
% L. Si, X. Zhang, Y. Zhang, Y. Tian, and S. Yang. Reinforcement learning-
% assisted multi-stage evolutionary constrained multi-objective
% optimization. ACM Transactions on Evolutionary Learning, 2025.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Initialization
            ArcSize = 2*Problem.N;
            ArcDec  = unifrnd(repmat(Problem.lower, ArcSize, 1), repmat(Problem.upper, ArcSize, 1)); 
            Archive = Problem.Evaluation(ArcDec);

            %% DQN model
            Learning   = 0.2;
            ModelExist = false;
            Counter    = 1;
            Greedy     = 0.1;
            Gama       = 0.9;
            Records    = [];
            ReSize     = 20;

            %% Action: 1-CDP, 2-MOP, 3-CDP+MOP
            ActionNum  = 2;
            StateNum   = 3;
            Iter       = 5;
            UpdateIter = 20;
            % Initial action
            Action = 1;

            %% Optimization
            while Algorithm.NotTerminated(Archive)
                %% Excute the action in iter iterations
                [ArchiveU, ratio] = ExcuteAction(Problem, Archive, Action, Iter);
             
                %% Recording
                % State of old Archive
                State = EstimateState(Archive,StateNum);

                % State of new Archive
                StateU = EstimateState(ArchiveU,StateNum);
                % Reward
                Reward = EstimateReward(Problem, Archive, ArchiveU, ratio,0,0,0);
                % Record
                Record = [State,Action,Reward,StateU];
                
                % Update Archive
                Archive = ArchiveU;

                % Update Records
                Records(mode(Counter-1,ReSize)+1,:) = Record;
                Counter = Counter + 1;

                %% Action selection
                if Problem.FE < Learning*Problem.maxFE
                    % Learning process
                    % Data collection
                    if Problem.FE/Problem.maxFE < Learning*0.2
                        % Random selection
                        Action = randi(ActionNum);
                    else
                        % DQN based selection
                        if ~ModelExist
                            % Build DQN model
                            [Model,Paras] = BuildRL(Records,StateNum);
                            ModelExist    = true;
                            Action        = randi(ActionNum);
                        else
                            if rand() < Greedy
                                % Random selection
                                Action = randi(ActionNum);
                            else
                                % DQN based selection
                                [Prob1,Prob2] = ApplyRL(State, Model, Paras);
                                Probs         = [Prob1;Prob2];
                                [~,Action]    = max(Probs(:,1));
                            end
                        end
                    end
                else
                    % Applying process
                    if rand() < Greedy
                        % Random selection
                        Action = randi(ActionNum);
                    else
                        % DQN based selection
                        [Prob1,Prob2] = ApplyRL(State, Model, Paras);
                        Probs         = [Prob1;Prob2];
                        [~,Action]    = max(Probs(:,1));
                    end
                end
                % Update model every UpdateIter iterations
                if ModelExist
                    if mode(Counter, UpdateIter) == 0
                        [Model, Paras] = UpdateRL(Records, Model, Paras,Gama, StateNum);
                    end
                end
            end
        end
    end
end