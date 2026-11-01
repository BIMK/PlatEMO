classdef SCMOEA < ALGORITHM
% <2026> <multi> <real/integer/binary> <large/none> <constrained/none> <sparse>
% Sparse constrained multi-objective evolutionary algorithm

%------------------------------- Reference --------------------------------
% S. Shao, Y. Tian, Y. Zhang, and X. Zhang. An evolutionary algorithm with
% adaptive search space switching and constraint-guided operators for
% sparse constrained multiobjective optimization. IEEE Transactions on
% Systems, Man, and Cybernetics: Systems, 2026.
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
            %% Population initialization
            % Calculate the fitness of each decision variable
            TDec    = [];
            TMask   = [];
            TempPop = [];
            Fitness = zeros(1,Problem.D);
            for i = 1 : 1+4*any(Problem.encoding~=4)
                Dec = unifrnd(repmat(Problem.lower,Problem.D,1),repmat(Problem.upper,Problem.D,1));
                Dec(:,Problem.encoding==4) = 1;
                Mask       = eye(Problem.D);
                Population = Problem.Evaluation(Dec.*Mask);
                TDec       = [TDec;Dec];
                TMask      = [TMask;Mask];
                TempPop    = [TempPop,Population];
                Fitness    = Fitness + NDSort([Population.objs,Population.cons],inf);
            end
            % Generate initial population
            Dec = unifrnd(repmat(Problem.lower,Problem.N,1),repmat(Problem.upper,Problem.N,1));
            Dec(:,Problem.encoding==4) = 1;
            Mask = false(Problem.N,Problem.D);
            for i = 1 : Problem.N
                Mask(i,TournamentSelection(2,ceil(rand*Problem.D),Fitness)) = 1;
            end
            Population  = Problem.Evaluation(Dec.*Mask);
            [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection([Population,TempPop],[Dec;TDec],[Mask;TMask],Problem.N,1);           
            Archive     = Population;
            LastArchive = Archive;
            [state, num_actions, ~, ~, ~] = GenerateSample(Problem, -1, LastArchive, Archive);
            num_states        = length(state);
            exp_replay_freq   = 1;
            copy_weights_freq = exp_replay_freq;
            batch_size        = 16;
            ddqn              = DDQNagent(num_states, num_actions, [16, 16], Problem, 1000);

            %% Optimization
            while Algorithm.NotTerminated(Archive)              
                action        = ddqn.action(state);                        
                ddqn.iter_num = ddqn.iter_num + 1;
                MatingPool    = TournamentSelection(2,2*Problem.N,FrontNo,-CrowdDis);
                if action == 1
                    [OffDec,OffMask] = OperatorSC1(Problem,Dec(MatingPool,:),Mask(MatingPool,:),Fitness,Mask);
                    Offspring        = Problem.Evaluation(OffDec.*OffMask);
                    [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Problem.N, 1);
                elseif action == 2
                    [OffDec,OffMask] = OperatorSC2(Problem,Dec(MatingPool,:),Mask(MatingPool,:),Fitness);
                    Offspring        = Problem.Evaluation(OffDec.*OffMask);
                    [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Problem.N, 1);
                elseif action == 3
                    [OffDec,OffMask] = OperatorSC3(Problem,Dec(MatingPool,:),Mask(MatingPool,:),Fitness);
                    Offspring        = Problem.Evaluation(OffDec.*OffMask);
                    [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Problem.N, 2);
                end
                Archive = NSGAII_EnvironmentalSelection([Archive,Population], Problem.N);
                [state, action, next_state, ~, reward] = GenerateSample(Problem, action, LastArchive, Archive);
                ddqn.store(state, action, reward, next_state);
                state       = next_state;
                LastArchive = Archive;
                if mod(ddqn.iter_num, exp_replay_freq) == 0
                    ddqn.experience_replay(batch_size);
                end
                if mod(ddqn.iter_num, copy_weights_freq) == 0
                    ddqn.copy_weights_agent_to_target(); 
                end
            end
        end
    end
end