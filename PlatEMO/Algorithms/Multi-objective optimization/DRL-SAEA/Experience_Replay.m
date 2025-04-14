classdef Experience_Replay < handle
    %EXPERIENCE_REPLAY 
    % A circular buffer
    properties
        buffer; % matrix array storing all of the experiences
        index = 1; % Current index of the latest insert
        overflow = 0; % Flag to determine whether the bugfer is full
        N; % Maximum memory limit
        nS; % Number of State 
    end
    
    methods
        % Constructor to initialise the buffer memory (is faster)
        function obj = Experience_Replay(n_S, maxER)
            obj.N = maxER;
            obj.buffer = zeros(maxER, (2*n_S + 2));
            obj.nS = n_S;
            
        end
        
        % Inserts the latest experience into to buffer 
        function insert_experience(obj,S,A,R,Snew)
            if obj.index > obj.N
                obj.index = 1;
                obj.overflow = 1;
            end
            nS = obj.nS;
            insert = zeros(1,(2*nS + 2));
            for i = 1:(2*nS + 2)
                if i < nS+1
                    insert(i) = S(i);
                elseif i ==  nS + 1 
                    insert(i) = A; 
                elseif i == nS + 2
                    insert(i) = R; 
                elseif i > (nS + 2)
                    insert(i) = Snew(i - nS - 2);
                end
            end
            obj.buffer(obj.index, 1:(2*nS + 2)) = insert;
            obj.index = obj.index+1;
        end
        
        % Samples a batch from the buffer with uniofrm random distribution
        function [Sold, A, R, Snew] = get_batch(obj, batch_size)
            if ~obj.overflow
                rand_index = obj.index-1;
            else 
                rand_index = obj.N;
            end
            all_i = randi(rand_index, [batch_size, 1]);
            batch = obj.buffer(all_i, 1:(2*obj.nS + 2));
            Sold = batch(1:batch_size, 1:obj.nS);
            A = batch(1:batch_size, obj.nS+1); 
            R = batch(1:batch_size, obj.nS+2);
            Snew = batch(1:batch_size, obj.nS+3:2*obj.nS+2); 
        end
    end
end


