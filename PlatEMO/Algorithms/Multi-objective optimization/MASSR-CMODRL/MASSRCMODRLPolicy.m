classdef MASSRCMODRLPolicy < handle
%MASSRCMODRLPolicy - Hybrid configuration policies for CMODRL.
%
%   Two agents use separate 16-dimensional MASSR states. A continuous actor
%   controls the normalized environment factor; a discrete Q-network chooses
%   DE or GA. Both learn from the same archive-improvement reward.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        Model;                      % Online actor, critic, and operator Q-network
        Target;                     % Target copies for bootstrap estimates
        Optimizer = [];             % Adam moments for the online parameters
        Buffer    = cell(1,200);    % Cyclic replay with 200 joint-transition slots
        Count     = 0;              % Number of stored joint transitions
        Position  = 0;              % Current replay write position
    end
    methods
        function obj = MASSRCMODRLPolicy()
        %MASSRCMODRLPolicy - Initialize environment and operator learners.
        %
        %   Create the continuous actor, its state-action critic, and the
        %   two-action operator Q-network. Initialize target networks from
        %   the online models.

            obj.Model.actor    = MASSRNN.initialize([16,30,30,1]);
            obj.Model.critic   = MASSRNN.initialize([17,30,30,1]);
            obj.Model.operator = MASSRNN.initialize([16,30,30,2]);
            obj.Target         = obj.Model;
        end
        function action = act(obj,state)
        %act - Select normalized environment and operator actions.
        %
        %   action(1) is in [0,1] and is mapped to [0,meanCV] by the main
        %   algorithm. action(2) selects DE or GA. Each learner receives its
        %   own 16-dimensional MASSR state.

            if obj.Count < 32 || rand < 0.05
                action = [rand;randi(2)];
            else
                u            = MASSRNN.forward(obj.Model.actor,state(:,1));
                q            = MASSRNN.forward(obj.Model.operator,state(:,2));
                [~,operator] = max(q);
                action       = [0.5*(1+tanh(u/2));operator];
            end
        end
        function remember(obj,state,action,reward,next,done)
        %remember - Store one shared-reward transition.
        %
        %   Store the two agent states, joint action, archive-improvement
        %   reward, next states, and terminal flag in a 200-slot cyclic
        %   buffer.

            obj.Position             = mod(obj.Position,200) + 1;
            obj.Count                = min(obj.Count+1,200);
            obj.Buffer{obj.Position} = struct('state',state,'action',action, ...
                'reward',reward,'next',next,'done',done);
        end
        function train(obj)
        %train - Update both learners from shared replay.
        %
        %   After 32 transitions, sample up to 60 records. Update the
        %   critic, deterministic continuous actor, and Double-Q operator
        %   learner; then softly update their target networks.

            if obj.Count < 32
                return;
            end

            % Sample joint transitions for both learners
            B     = min(60,obj.Count);
            batch = obj.Buffer(randperm(obj.Count,B));
            S     = zeros(16,2,B);
            T     = S;
            A     = zeros(2,B);
            R     = zeros(1,B);
            done  = R;
            for i = 1 : B
                S(:,:,i) = batch{i}.state;
                T(:,:,i) = batch{i}.next;
                A(:,i)   = batch{i}.action;
                R(i)     = batch{i}.reward;
                done(i)  = batch{i}.done;
            end
            s1 = reshape(S(:,1,:),16,B);
            s2 = reshape(S(:,2,:),16,B);
            t1 = reshape(T(:,1,:),16,B);
            t2 = reshape(T(:,2,:),16,B);

            % Update the continuous environment-factor critic and actor
            nextU           = 0.5*(1+tanh(MASSRNN.forward(obj.Target.actor,t1)/2));
            nextQ           = MASSRNN.forward(obj.Target.critic,[t1;nextU]);
            y               = R + 0.95*(1-done).*nextQ;
            [q,c]           = MASSRNN.forward(obj.Model.critic,[s1;A(1,:)]);
            G.critic        = MASSRNN.backward(obj.Model.critic,c,2*(q-y)/B);
            [u,c]           = MASSRNN.forward(obj.Model.actor,s1);
            u               = 0.5*(1+tanh(u/2));
            [~,criticCache] = MASSRNN.forward(obj.Model.critic,[s1;u]);
            [~,dInput]      = MASSRNN.backward(obj.Model.critic,criticCache,-ones(1,B)/B);
            G.actor         = MASSRNN.backward(obj.Model.actor,c,dInput(end,:).*u.*(1-u));

            % Update the discrete DE/GA learner with a Double-Q target
            [q,c]      = MASSRNN.forward(obj.Model.operator,s2);
            online     = MASSRNN.forward(obj.Model.operator,t2);
            future     = MASSRNN.forward(obj.Target.operator,t2);
            [~,best]   = max(online,[],1);
            nextQ      = future(sub2ind(size(future),best,1:B));
            y          = R + 0.95*(1-done).*nextQ;
            index      = sub2ind(size(q),A(2,:),1:B);
            dq         = zeros(size(q));
            dq(index)  = 2*(q(index)-y)/B;
            G.operator = MASSRNN.backward(obj.Model.operator,c,dq);

            % Apply the joint update and track target networks
            [obj.Model,obj.Optimizer] = MASSRNN.adam(obj.Model,G,obj.Optimizer,1e-3);
            obj.Target                = MASSRNN.add(MASSRNN.scale(obj.Target,0.99),MASSRNN.scale(obj.Model,0.01));
        end
    end
end