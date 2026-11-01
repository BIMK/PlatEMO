classdef MASSRMADEPolicy < handle
%MASSRMADEPolicy - Value-decomposition policies for MADE.
%
%   Three agents select the mating strategy, DE operator, and F/CR setting.
%   Their independent Q-networks share one team temporal-difference
%   objective and use fixed MASSR states.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        Model;          % Online Q-networks for the three agents
        Dimensions;     % Numbers of available actions per agent
        Target;         % Target Q-networks for bootstrapping
        Updates = 0;    % Number of completed gradient updates
        Count   = 0;    % Number of stored joint transitions
    end
    properties(Access = private)
        Buffer;                   % Cyclic replay of joint transitions
        Position       = 0;       % Current replay write position
        Optimizer      = [];      % Adam moments for the online networks
        RewardMean     = 0;       % Running mean of team rewards
        RewardVariance = 1;       % Running variance of team rewards
        RewardCount    = 1e-4;    % Sample count for running reward moments
    end

    methods
        function obj = MASSRMADEPolicy()
        %MASSRMADEPolicy - Initialize the agent Q-networks and replay buffer.
        %
        %   Create 3 independent 16-input Q-networks, copy their target
        %   networks, and reserve 5000 joint-transition replay slots.

            dimensions     = [3,4,4];
            obj.Dimensions = dimensions;
            obj.Model      = cell(1,numel(dimensions));
            for k = 1 : numel(obj.Dimensions)
                obj.Model{k} = MASSRNN.initialize([16,64,64,dimensions(k)]);
            end
            obj.Target = obj.Model;
            obj.Buffer = cell(1,5000);
        end

        function action = act(obj,state,available,progress)
        %act - Select a joint action with epsilon-greedy exploration.
        %
        %   state is 16-by-3; available marks valid discrete actions.
        %   progress is the normalized policy-learning progress used to
        %   reduce epsilon from 1 to 0.05.

            epsilon = 1 - 0.95*min(max(progress,0),1);
            action  = zeros(numel(obj.Dimensions),1);
            for k = 1 : numel(obj.Dimensions)
                q     = MASSRNN.forward(obj.Model{k},state(:,k));
                valid = find(available(1:numel(q),k));
                if rand < epsilon
                    action(k) = valid(randi(numel(valid)));
                else
                    [~,best]  = max(q(valid));
                    action(k) = valid(best);
                end
            end
        end

        function remember(obj,state,action,reward,nextState,done,nextAvailable)
        %remember - Append one joint transition to cyclic replay.
        %
        %   Store all agent states and actions, the scalar team reward, the
        %   next states, the terminal flag, and the next-action availability
        %   mask.

            obj.Position             = mod(obj.Position,5000) + 1;
            obj.Count                = min(obj.Count+1,5000);
            obj.Buffer{obj.Position} = struct('state',state,'action',action(:), ...
                'reward',reward,'nextState',nextState,'done',logical(done),'available',nextAvailable);
        end

        function train(obj)
        %train - Perform one minibatch update when replay is ready.
        %
        %   Sample 32 joint transitions, standardize rewards, and optimize
        %   the summed team Q-value. Target networks synchronize every 200
        %   updates.

            if obj.Count < 32
                return;
            end

            % Assemble a minibatch of joint transitions
            records         = obj.Buffer(randperm(obj.Count,32));
            values          = recordsToCell(records,'state');
            batch.state     = cat(3,values{:});
            values          = recordsToCell(records,'nextState');
            batch.nextState = cat(3,values{:});
            values          = recordsToCell(records,'action');
            batch.action    = cat(2,values{:});
            values          = recordsToCell(records,'reward');
            batch.reward    = cat(2,values{:});
            values          = recordsToCell(records,'done');
            batch.done      = cat(2,values{:});
            values          = recordsToCell(records,'available');
            batch.available = cat(3,values{:});

            % Standardize the shared reward and update all local Q-networks
            obj.standardize(batch.reward);
            batch.reward              = (batch.reward-obj.RewardMean)/sqrt(obj.RewardVariance+1e-8);
            [~,gradient]              = MASSRMADEPolicy.objective(obj.Model,obj.Target,batch);
            [obj.Model,obj.Optimizer] = MASSRNN.adam(obj.Model,gradient,obj.Optimizer,1e-4,10);
            obj.Updates               = obj.Updates + 1;
            if mod(obj.Updates,200) == 0
                obj.Target = obj.Model;
            end
        end
    end

    methods(Access = private)
        function standardize(obj,reward)
        %standardize - Update running team-reward moments.
        %
        %   Combine the minibatch mean and population variance with the
        %   accumulated running statistics.

            n                  = numel(reward);
            delta              = mean(reward) - obj.RewardMean;
            total              = obj.RewardCount + n;
            varianceSum        = obj.RewardVariance*obj.RewardCount + var(reward,1)*n + ...
                delta^2*obj.RewardCount*n/total;
            obj.RewardMean     = obj.RewardMean + delta*n/total;
            obj.RewardVariance = varianceSum/total;
            obj.RewardCount    = total;
        end
    end

    methods(Static)
        function [loss,G] = objective(P,target,batch)
        %objective - Evaluate the VDN team TD loss and gradients.
        %
        %   The joint Q-value is the sum of selected local Q-values. Online
        %   networks select available next actions; target networks evaluate
        %   them. Every selected local value receives the shared team TD
        %   derivative.

            B         = numel(batch.reward);
            chosen    = zeros(1,B);
            bootstrap = zeros(1,B);
            cache     = cell(size(P));
            q         = cell(size(P));
            indices   = cell(size(P));

            % Sum selected local values and masked bootstrap estimates
            for k = 1 : numel(P)
                state           = reshape(batch.state(:,k,:),16,B);
                next            = reshape(batch.nextState(:,k,:),16,B);
                [q{k},cache{k}] = MASSRNN.forward(P{k},state);
                rows            = batch.action(k,:);
                indices{k}      = sub2ind(size(q{k}),rows,1:B);
                chosen          = chosen + q{k}(indices{k});
                online          = MASSRNN.forward(P{k},next);
                future          = MASSRNN.forward(target{k},next);
                for i = 1 : B
                    if batch.done(i)
                        continue;
                    end
                    valid        = find(batch.available(1:size(online,1),k,i));
                    [~,best]     = max(online(valid,i));
                    bootstrap(i) = bootstrap(i) + future(valid(best),i);
                end
            end

            % The shared team TD error supervises every selected local value
            expected = batch.reward + 0.99*(1-double(batch.done)).*bootstrap;
            delta    = chosen - expected;
            loss     = mean(delta.^2);
            G        = cell(size(P));
            for k = 1 : numel(P)
                dQ             = zeros(size(q{k}));
                dQ(indices{k}) = 2*delta/B;
                G{k}           = MASSRNN.backward(P{k},cache{k},dQ);
            end
        end
    end
end

function values = recordsToCell(records,key)
%recordsToCell - Extract one field from replay records.

    values = cell(size(records));
    for i = 1 : numel(records)
        values{i} = records{i}.(key);
    end
end