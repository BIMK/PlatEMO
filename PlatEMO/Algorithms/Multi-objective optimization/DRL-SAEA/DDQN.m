classdef DDQN < handle   
    properties
        number_state;
        number_actions;
        memory;     
        epsilon_decay_start; 
        epsilon_decay_finish; 
        epsilon_decay_rate;
        discount_factor = 0.999; 
        agent;
        target;
        loss = 0;
        reward_W = [1,1,1,1]
        epoch = 0;
        max_epochs = 10000;
        training_records;
        action_count = 0;     
    end
    
    methods
        function obj = DDQN(n_S, n_A, layers, maxER)           
            obj.number_state = n_S;
            obj.number_actions = n_A;
            obj.memory = Experience_Replay(n_S, maxER);
            obj.agent = makenet(layers);
            obj.target = makenet(layers);
            obj.copy_weights_agent_to_target()
            obj.agent = confignet(obj.agent, rand(n_S, 5), -rand(n_A, 5));
            obj.target = confignet(obj.target, rand(n_S, 5), -rand(n_A, 5));
            obj.agent.trainParam.showWindow=0;
            obj.target.trainParam.showWindow=0;
            obj.training_records = zeros(0, obj.max_epochs);            
        end
                   
        function copy_weights_agent_to_target(obj)
            obj.target.IW = obj.agent.IW; 
            obj.target.LW = obj.agent.LW;
            obj.target.b = obj.agent.b;
        end
                
        function a = action(obj, state)            
            action_vals = obj.agent(transpose(state));           
            a = datasample([repmat(find(action_vals==max(action_vals)),1,obj.number_actions) 1:obj.number_actions],1);            
        end
                
        function a = exploit_action(obj, state)
            action_vals = obj.agent(transpose(state));
            [~, a] = max(action_vals);
        end
        
        function store(obj, S, A, R, S1)
            obj.memory.insert_experience(S,A,R,S1);
        end
                  
        function experience_replay(obj, batch_size) 
           obj.epoch = obj.epoch+1;
           [S, A, R, Sn] = obj.memory.get_batch(batch_size);  
           st_predict = obj.agent(transpose(S));
           nst_predict = obj.agent(transpose(Sn));
           nst_predict_target = obj.target(transpose(Sn));
           Q_target = st_predict;           
           for i = 1:batch_size
               [~,next_best_action] = max(nst_predict(:,i));
               Q_target(A(i),i) = R(i) + obj.discount_factor * nst_predict_target(next_best_action,i);
           end
           obj.agent = trainnet(obj.agent, transpose(S), Q_target);        
           obj.epoch = obj.epoch + 1;
        end
    end
end


function net = makenet(layers)
    net = feedforwardnet(layers);
    net.layers{1:end-1}.transferFcn = 'tansig'; 
    net.layers{end}.transferFcn = 'purelin';
end

function net = confignet(net, x, t)
    net = configure(net, x, t);
end

function net = trainnet(net, x, t)
        net = train(net, x, t, 'useParallel','no');
end
