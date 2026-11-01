classdef MASSRMARLACSPolicy < handle
%MASSRMARLACSPolicy - Three hybrid actors with a centralized critic.
%
%   Three actors control the guiding criterion, normalized selection ratio,
%   and sparse operator. Each receives one 16-dimensional MASSR state. The
%   critic receives all 48 state values and the six-component joint action.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        Model;                       % Three online actors and a centralized critic
        Target;                      % Target actors and centralized critic
        Optimizer = [];              % Adam moments for the online parameters
        Buffer    = cell(1,1000);    % Cyclic replay with 1000 joint-transition slots
        Count     = 0;               % Number of stored joint transitions
        Position  = 0;               % Current replay write position
    end
    methods
        function obj = MASSRMARLACSPolicy()
        %MASSRMARLACSPolicy - Initialize three actors and the joint critic.
        %
        %   The actors output three guiding logits, one normalized ratio,
        %   and two operator logits. The critic takes 48 state values plus
        %   six joint-action values.

            widths = [3,1,2];
            for k = 1 : 3
                obj.Model.actor{k} = MASSRNN.initialize([16,64,32,widths(k)]);
            end
            obj.Model.critic = MASSRNN.initialize([54,256,128,1]);
            obj.Target       = obj.Model;
        end
        function [action,vector] = act(obj,state,progress)
        %act - Select a hybrid action and its critic encoding.
        %
        %   state is 16-by-3. action contains a guiding index, a ratio in
        %   [0,1], and an operator index. vector concatenates three guiding
        %   indicators, the ratio, and two operator indicators.

            if obj.Count < 64 || rand < 0.3^progress
                action = [randi(3);rand;randi(2)];
            else
                [vector,~] = actorOutput(obj.Model.actor,state(:));
                [~,a]      = max(vector(1:3));
                [~,b]      = max(vector(5:6));
                action     = [a;vector(4);b];
            end
            vector              = zeros(6,1);
            vector(action(1))   = 1;
            vector(4)           = action(2);
            vector(4+action(3)) = 1;
        end
        function remember(obj,state,vector,reward,next,done)
        %remember - Store a joint transition for centralized training.
        %
        %   Flatten the three states into 48 values and store the six-value
        %   action encoding with the scalar HV reward, next states, and
        %   terminal flag.

            obj.Position             = mod(obj.Position,1000) + 1;
            obj.Count                = min(obj.Count+1,1000);
            obj.Buffer{obj.Position} = struct('state',state(:),'action',vector, ...
                'reward',reward,'next',next(:),'done',done);
        end
        function train(obj)
        %train - Update the centralized critic and three actors.
        %
        %   Sample 64 joint transitions. Backpropagate the critic action
        %   gradient through the continuous sigmoid and the discrete softmax
        %   outputs, then softly update all target networks.

            if obj.Count < 64
                return;
            end

            % Sample centralized state and action transitions
            B       = 64;
            records = obj.Buffer(randperm(obj.Count,B));
            S       = zeros(48,B);
            T       = S;
            A       = zeros(6,B);
            R       = zeros(1,B);
            done    = R;
            for i = 1 : B
                S(:,i)  = records{i}.state;
                T(:,i)  = records{i}.next;
                A(:,i)  = records{i}.action;
                R(i)    = records{i}.reward;
                done(i) = records{i}.done;
            end

            % Regress the critic to the target-policy return
            next          = actorOutput(obj.Target.actor,T);
            future        = MASSRNN.forward(obj.Target.critic,[T;next]);
            target        = R + 0.99*(1-done).*future;
            [value,cache] = MASSRNN.forward(obj.Model.critic,[S;A]);
            G.critic      = MASSRNN.backward(obj.Model.critic,cache,2*(value-target)/B);

            % Backpropagate the critic action gradient through each actor
            [current,cache] = actorOutput(obj.Model.actor,S);
            [~,c]           = MASSRNN.forward(obj.Model.critic,[S;current]);
            [~,dInput]      = MASSRNN.backward(obj.Model.critic,c,-ones(1,B)/B);
            dA              = dInput(49:end,:);
            ranges          = {1:3,4,5:6};
            for k = 1 : 3
                rows       = ranges{k};
                output     = current(rows,:);
                derivative = dA(rows,:);
                if k == 2
                    derivative = derivative.*output.*(1-output);
                else
                    derivative = output.*(derivative-sum(derivative.*output,1));
                end
                G.actor{k} = MASSRNN.backward(obj.Model.actor{k},cache{k},derivative);
            end

            % Apply the joint update and softly track target networks
            [obj.Model,obj.Optimizer] = MASSRNN.adam(obj.Model,G,obj.Optimizer,1e-3);
            obj.Target                = MASSRNN.add(MASSRNN.scale(obj.Target,0.99),MASSRNN.scale(obj.Model,0.01));
        end
    end
end

function [actions,cache] = actorOutput(actors,S)
%actorOutput - Evaluate the three differentiable actor outputs.
%
%   S contains 48 stacked state values per column. Concatenate a three-way
%   softmax, a sigmoid ratio, and a two-way softmax into six action values.

    B       = size(S,2);
    actions = zeros(6,B);
    cache   = cell(1,3);
    ranges  = {1:3,4,5:6};
    for k = 1 : 3
        [output,cache{k}] = MASSRNN.forward(actors{k},S((k-1)*16+(1:16),:));
        if k == 2
            output = 0.5*(1+tanh(output/2));
        else
            output = exp(output-max(output,[],1));
            output = output./sum(output,1);
        end
        actions(ranges{k},:) = output;
    end
end