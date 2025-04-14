classdef DRLSAEA < ALGORITHM
% <2025> <multi> <real> <expensive> <constrained>
% Deep reinforcement learning-based expensive constrained evolutionary algorithm
% wmax --- 20 --- Number of generations before updating surrogate models
% mu   ---  5 --- Number of real evaluated solutions at each iteration

%------------------------------- Reference --------------------------------
% S. Shao, Y. Tian, and Y. Zhang. Deep reinforcement learning assisted
% surrogate model management for expensive constrained multi-objective
% optimization. Swarm and Evolutionary Computation, 2025, 92: 101817.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)           
             %% Parameter setting
            [wmax,mu] = Algorithm.ParameterSet(20,5);
            
            %% Generate the initial population 
            NI          = 11*Problem.D-1;
            P           = UniformPoint(NI,Problem.D,'Latin');
            Population  = Problem.Evaluation(repmat(Problem.upper-Problem.lower,NI,1).*P+repmat(Problem.lower,NI,1));                                            
            Archive     = UpdateArchive(Population,Problem.N);
            THETA_OBJ   = 5.*ones(Problem.M,Problem.D);
            THETA_CV    = 5.*ones(1,Problem.D);
            THETA_CON   = 5.*ones(size(Population.cons,2),Problem.D);
            
            LastArchive = Archive;
            [state, num_actions, ~, ~, ~] = GenerateSample(Problem, -1, LastArchive, Archive);
            num_states = length(state);
            exp_replay_freq = 8;
            copy_weights_freq = exp_replay_freq + 3;
            batch_size = 16;
            step = 0;
            ddqn = DDQN(num_states, num_actions, [5, 5], 100);

            %% Optimization
            while Algorithm.NotTerminated(Archive) 
                action = ddqn.action(state);           
                step = step + 1;
                if action == 1
                    status  = 1;
                    CV      = NormalizeCV(Population.cons); 
                    PopDec  = Population.decs;
                    PopObj  = [Population.objs,CV];
                    M       = Problem.M+1;
                    Model   = cell(1,M);
                    Fitness = CalFitness(Population.objs,Population.cons);
                    THETA   = [THETA_OBJ;THETA_CV];              
                    for i = 1 : M
                        dmodel     = dacefit(PopDec,PopObj(:,i),'regpoly0','corrgauss',THETA(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                        Model{i}   = dmodel;
                        THETA(i,:) = dmodel.theta;
                    end
                    THETA_OBJ = THETA(1:Problem.M,:);
                    THETA_CV  = THETA(end,:);
                elseif action == 2
                    status   = 2;
                    PopCon   = max(0,Population.cons);
                    MaxCon   = max(PopCon);
                    index    = find(MaxCon>0);
                    PopCon(:,MaxCon==0) = [];
                    PopCon   = (PopCon-min(PopCon,[],1))./(max(PopCon,[],1)-min(PopCon,[],1));
                    PopDec   = Population.decs;
                    PopObj   = [Population.objs,PopCon];
                    M        = size(PopObj,2);
                    Model    = cell(1,M);
                    Fitness  = CalFitness([Population.objs,sum(max(0,Population.cons),2)]);
                    THETA    = [THETA_OBJ;THETA_CON(index,:)];                  
                    for i = 1 : M
                        dmodel     = dacefit(PopDec,PopObj(:,i),'regpoly0','corrgauss',THETA(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                        Model{i}   = dmodel;
                        THETA(i,:) = dmodel.theta;
                    end
                    THETA_OBJ = THETA(1:Problem.M,:);
                    THETA_CON(index,:) = THETA(Problem.M+1:end,:);
                elseif action == 3
                    status   = 3;
                    PopDec   = Population.decs;
                    PopObj   = Population.objs;
                    M        = Problem.M;
                    Model    = cell(1,M);
                    Fitness  = CalFitness(PopObj);                    
                    for i = 1 : M
                        dmodel     = dacefit(PopDec,PopObj(:,i),'regpoly0','corrgauss',THETA_OBJ(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                        Model{i}   = dmodel;
                        THETA_OBJ(i,:) = dmodel.theta;
                    end
                end             
                for w = 1: wmax
                    drawnow();
                    MatingPool = TournamentSelection(2,NI,Fitness);
                    OffDec     = OperatorGA(Problem,PopDec(MatingPool,:));
                    PopDec     = cat(1,PopDec,OffDec);
                    [N,~]      = size(PopDec);
                    PopObj     = zeros(N,M);
                    MSE        = zeros(N,M);
                    for i = 1: N
                        for j = 1 : M
                            [PopObj(i,j),~,MSE(i,j)] = predictor(PopDec(i,:),Model{j});
                        end
                    end
                    [PopDec,PopObj,Fitness] = EnvironmentalSelection(PopDec,PopObj,NI,Problem.M,status);
                end               
                [NewDec,~,~] = EnvironmentalSelection(PopDec,PopObj,mu,Problem.M,status);
                New          = Problem.Evaluation(NewDec);              
                Population   = UpdatePopulation(Population,New,NI-mu,status);
                Archive      = cat(2,Archive,New);
                Archive      = UpdateArchive(Archive,Problem.N);
                [state, action, next_state, ~, reward] = GenerateSample(Problem, action, LastArchive, Archive);
                ddqn.store(state, action, reward, next_state);
                state = next_state;
                LastArchive = Archive;               
                if mod(step, exp_replay_freq) == 0
                    ddqn.experience_replay(batch_size);
                end
                if mod(step, copy_weights_freq) == 0
                    ddqn.copy_weights_agent_to_target(); 
                end
            end          
        end
    end  
end

function CV = NormalizeCV(PopCon)
    PopCon = max(0,PopCon);
    PopCon = (PopCon-min(PopCon,[],1))./(max(PopCon,[],1)-min(PopCon,[],1));
    PopCon(:,isnan(PopCon(1,:))) = 0;
    CV = sum(PopCon,2);
end