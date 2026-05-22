classdef DDPG < handle

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        state_dim;
        action_dim;
        action_bound;
        opts;
        
        actor;
        actor_target;
        critic;
        critic_target;
        
        ou_theta;
        ou_sigma;
        ou_mu;
        ou_state;

        avgGradCritic;
        avgSqCritic;
        avgGradActor;
        avgSqActor;
        adamIter;

        buffer_capacity;
        buffer_ptr;
        buffer_size_count;
        buf_states;
        buf_actions;
        buf_rewards;
        buf_next_states;
        buf_dones;
    end
    
    methods
        function obj = DDPG(state_dim, action_dim, action_bound, optOverrides)
            if nargin < 4
                optOverrides = struct();
            end
            obj.state_dim    = state_dim;
            obj.action_dim   = action_dim;
            obj.action_bound = action_bound;
            obj.opts         = DDPG.mergeOptions(DDPG.defaultOptions(), optOverrides);
            
            buffer_cap   = obj.opts.buffer_size;
            hidden_sizes = obj.opts.hidden_sizes;
            
            obj.ou_state = zeros(1, action_dim);
            
            obj.avgGradCritic = [];
            obj.avgSqCritic   = [];
            obj.avgGradActor  = [];
            obj.avgSqActor    = [];
            obj.adamIter      = 0;
            
            obj.buffer_capacity   = buffer_cap;
            obj.buffer_ptr        = 1;
            obj.buffer_size_count = 0;
            obj.buf_states        = zeros(buffer_cap, state_dim);
            obj.buf_actions       = zeros(buffer_cap, action_dim);
            obj.buf_rewards       = zeros(buffer_cap, 1);
            obj.buf_next_states   = zeros(buffer_cap, state_dim);
            obj.buf_dones         = zeros(buffer_cap, 1);
            
            obj.actor = DDPG.buildActorDlnetwork(state_dim, action_dim, hidden_sizes);
            obj.actor_target = DDPG.buildActorDlnetwork(state_dim, action_dim, hidden_sizes);
            DDPG.copyLearnables(obj.actor_target, obj.actor);
            
            obj.critic = DDPG.buildCriticDlnetwork(state_dim, action_dim, hidden_sizes);
            obj.critic_target = DDPG.buildCriticDlnetwork(state_dim, action_dim, hidden_sizes);
            DDPG.copyLearnables(obj.critic_target, obj.critic);
        end
        
        function action = selectAction(obj, state, explore)
            if nargin < 3
                explore = true;
            end
            state  = single(reshape(state, [], 1));
            s_dl   = dlarray(state, 'CB');
            a_raw  = forward(obj.actor, s_dl);
            action = extractdata(a_raw)' * obj.action_bound;
            if explore
                ot = obj.opts.ou_theta;
                os = obj.opts.ou_sigma;
                om = obj.opts.ou_mu;
                obj.ou_state = obj.ou_state + ot * (om - obj.ou_state) + os * randn(1, obj.action_dim);
                noise  = obj.ou_state * obj.opts.noise_scale;
                action = action + noise;
                ab     = obj.action_bound;
                action = max(-ab, min(ab, action));
            end
        end
        
        function store(obj, state, action, reward, next_state, done)
            state      = reshape(state, 1, []);
            action     = reshape(action, 1, []);
            next_state = reshape(next_state, 1, []);
            p = obj.buffer_ptr;
            obj.buf_states(p, :)  = state;
            obj.buf_actions(p, :) = action;
            obj.buf_rewards(p)    = reward;
            obj.buf_next_states(p, :) = next_state;
            obj.buf_dones(p) = done;
            obj.buffer_ptr   = mod(p, obj.buffer_capacity) + 1;
            obj.buffer_size_count = min(obj.buffer_size_count + 1, obj.buffer_capacity);
        end
        
        function [critic_loss, actor_loss] = train(obj)
            critic_loss = 0;
            actor_loss  = 0;
            bs          = obj.opts.batch_size;
            if obj.buffer_size_count < bs
                return;
            end
            n           = min(bs, obj.buffer_size_count);
            idx         = randi(obj.buffer_size_count, [n, 1]);
            states      = obj.buf_states(idx, :);
            actions     = obj.buf_actions(idx, :);
            rewards     = obj.buf_rewards(idx);
            next_states = obj.buf_next_states(idx, :);
            dones       = obj.buf_dones(idx);
            
            s_next      = single(next_states');
            s_dl_next   = dlarray(s_next, 'CB');
            a_next_tanh = forward(obj.actor_target, s_dl_next);
            a_next      = extractdata(a_next_tanh) * obj.action_bound;
            sa_next     = [s_next; a_next];
            q_tgt_num   = extractdata(forward(obj.critic_target, dlarray(sa_next, 'CB')));
            y           = single(rewards + obj.opts.gamma .* (1 - dones) .* q_tgt_num(:));
            
            sa    = single([states'; actions']);
            y_row = reshape(y, 1, []);
            obj.adamIter = obj.adamIter + 1;
            [critic_loss, gradC] = dlfeval(@ddpgCriticGradients, obj.critic, dlarray(sa, 'CB'), dlarray(y_row, 'CB'));
            [obj.critic, obj.avgGradCritic, obj.avgSqCritic] = adamupdate(obj.critic, gradC, obj.avgGradCritic, obj.avgSqCritic, obj.adamIter, obj.opts.lr_critic);
            
            s_tr = single(states');
            [actor_loss, gradA] = dlfeval(@ddpgActorGradients, obj.actor, obj.critic, dlarray(s_tr, 'CB'), obj.action_bound);
            [obj.actor, obj.avgGradActor, obj.avgSqActor] = adamupdate(obj.actor, gradA, obj.avgGradActor, obj.avgSqActor, obj.adamIter, obj.opts.lr_actor);
            
            DDPG.softUpdateLearnables(obj.actor_target, obj.actor, obj.opts.tau);
            DDPG.softUpdateLearnables(obj.critic_target, obj.critic, obj.opts.tau);
            
            o             = obj.opts;
            o.noise_scale = o.noise_scale * o.noise_decay;
            obj.opts      = o;
        end
        
        function resetNoise(obj)
            obj.ou_state = zeros(1, obj.action_dim);
        end
    end
    
    methods (Static)
        function o = defaultOptions()         
            o              = struct();
            o.gamma        = 0.99;
            o.tau          = 0.005;
            o.lr_actor     = 1e-4;
            o.lr_critic    = 1e-3;
            o.batch_size   = 32;
            o.buffer_size  = 500;
            o.hidden_sizes = [32, 32];
            o.noise_scale  = 0.3;
            o.noise_decay  = 0.9995;
            o.ou_theta     = 0.15;
            o.ou_sigma     = 0.2;
            o.ou_mu        = 0;
        end
    end
    
    methods (Static, Access = private)
        function out = mergeOptions(base, override)
            out = base;
            if nargin < 2 || isempty(override)
                return;
            end
            fn = fieldnames(override);
            for i = 1 : numel(fn)
                f = fn{i};
                if isfield(base, f)
                    out.(f) = override.(f);
                end
            end
        end
        function net = buildActorDlnetwork(state_dim, action_dim, hidden_sizes)
            layers = featureInputLayer(state_dim, 'Normalization', 'none', 'Name', 'in');
            for k = 1 : numel(hidden_sizes)
                layers = [layers
                    fullyConnectedLayer(hidden_sizes(k), 'Name', sprintf('a_fc%d', k))
                    reluLayer('Name', sprintf('a_relu%d', k))];
            end
            layers = [layers
                fullyConnectedLayer(action_dim, 'Name', 'a_out_fc')
                tanhLayer('Name', 'a_out_tanh')];
            lg  = layerGraph(layers);
            net = dlnetwork(lg);
        end
        
        function net = buildCriticDlnetwork(state_dim, action_dim, hidden_sizes)
            in_dim = state_dim + action_dim;
            layers = featureInputLayer(in_dim, 'Normalization', 'none', 'Name', 'cin');
            for k = 1 : numel(hidden_sizes)
                layers = [layers
                    fullyConnectedLayer(hidden_sizes(k), 'Name', sprintf('c_fc%d', k))
                    reluLayer('Name', sprintf('c_relu%d', k))];
            end
            layers = [layers
                fullyConnectedLayer(1, 'Name', 'c_q')];
            lg  = layerGraph(layers);
            net = dlnetwork(lg);
        end
        
        function copyLearnables(dstNet, srcNet)
            Ld = dstNet.Learnables;
            Ls = srcNet.Learnables;
            for i = 1 : height(Ld)
                Ld.Value{i} = dlarray(extractdata(Ls.Value{i}));
            end
            dstNet.Learnables = Ld;
        end
        
        function softUpdateLearnables(targetNet, onlineNet, tau)
            Lt = targetNet.Learnables;
            Lo = onlineNet.Learnables;
            for i = 1 : height(Lt)
                Lt.Value{i} = tau * Lo.Value{i} + (1 - tau) * Lt.Value{i};
            end
            targetNet.Learnables = Lt;
        end
    end
end

function [loss, gradients] = ddpgCriticGradients(critic, sa_dl, y_dl)
    q         = forward(critic, sa_dl);
    loss      = mean((q - y_dl).^2, 'all');
    gradients = dlgradient(loss, critic.Learnables);
end

function [loss, gradients] = ddpgActorGradients(actor, critic, s_dl, action_bound)
    a_tanh    = forward(actor, s_dl);
    a         = a_tanh * action_bound;
    sa        = cat(1, s_dl, a);
    q         = forward(critic, sa);
    loss      = -mean(q, 'all');
    gradients = dlgradient(loss, actor.Learnables);
end