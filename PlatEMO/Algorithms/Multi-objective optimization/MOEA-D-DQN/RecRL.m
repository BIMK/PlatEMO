classdef RecRL < handle
% Choosing crossover operator via reinforcement learning

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        name              % operator name 
        problem         
        lambda_           % weight vector
        Opers             % candidate operators
        n                 % candidate operators number
        dqn               % Deep Q-Network 
        SW	              % sliding window
        a                 % action
        state             % current state
        state_            % next state
        countOpers        % count the selection frequency of different operators
        gen               % current genetic generation
    end
    methods
        function obj = RecRL(problem, lambda_, maxgen, NIND)
            obj.name       = 'RecRL';
            obj.problem    = problem;
            obj.lambda_    = lambda_;
            obj.Opers      = {Recsbx(), RecM2m(maxgen), DE_rand_1(), DE_rand_2()};
            obj.n          = length(obj.Opers);
            obj.dqn        = DQN(problem.D + problem.M, obj.n);
            obj.SW         = zeros(2, NIND * 4);
            obj.a          = 0;
            obj.state      = [];
            obj.state_     = [];
            obj.countOpers = zeros(obj.n, 1);
        end
        function [obj, offChrom] = do(obj, OldChrom, r0, neighbourVector, currentGen)
      	% Choose operator and generate offspring
            obj.gen   = currentGen;
            obj.state = [OldChrom(r0, :), obj.lambda_(r0, :)];
            if obj.dqn.memory_counter > 300
                obj.a = obj.dqn.choose_action(obj.state');
                for i = 1 : obj.n
                    if sum(obj.SW(1, :) == i) == 0
                        obj.a = i;
                        break;
                    end
                end
            else
                obj.a = randi(obj.n);
            end
            obj.countOpers(obj.a) = obj.countOpers(obj.a) + 1;
            if obj.a == 1
                offChrom = obj.Opers{1}.do(OldChrom, r0, neighbourVector);
            elseif obj.a == 2
                offChrom = obj.Opers{2}.do(OldChrom, r0, neighbourVector, currentGen);
            else
                offChrom = obj.Opers{obj.a}.do(OldChrom, r0, neighbourVector);
            end
            obj.state_ = [offChrom(1, :), obj.lambda_(r0, :)];
        end
        function obj = learn(obj, r)
        % Store state transition and (if necessary) update parameters in DQN
            obj.SW = [obj.SW(:, 2:end), [obj.a; r]];
            reward = max(obj.SW(2, obj.SW(1, :) == obj.a));
            obj.dqn.store_transition(obj.state, obj.a, reward, obj.state_);
            if obj.dqn.memory_counter > 200 && rand() < 0.2
                obj.dqn.learn();
            end
        end
    end
end