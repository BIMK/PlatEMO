 function Output = FCNForward(Input, instance, s_list)
% Fully connected neural network forward propagation
% The function input includes a set of neural network weights, 
% problem information, and a list of neural network structures
    
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    % Initialize output
    Output = zeros(size(Input, 1), size(instance, 1));
    for n = 1 : size(Input, 1)
        output = instance;
        % Initialize pointer
        pointer = 1;
        for i = 1 : size(s_list, 1)
            % Traverse the list of neural network structures
            if s_list(i, 2) ~= -1
                % Check if there is only bias, if not then
                weight  = Input(n, pointer:pointer + s_list(i, 1)*s_list(i, 2)-1);
                pointer = pointer + s_list(i, 1)*s_list(i, 2);
                output  = output * reshape(weight, [s_list(i, 1), s_list(i, 2)]);
            else
                % only bias then
                bias    = Input(n, pointer:pointer + s_list(i, 1)-1);
                pointer = pointer + s_list(i, 1);
                output  = output + bias;
                % Check if it is the last layer of the neural network
                if i == size(s_list, 1)
                    % If it is the last layer, output after activation
                    output       = step(output);
                    Output(n, :) = output;
                else
                    output = leaky_relu(output);
                end
            end
        end
    end
end

function y = leaky_relu(x)
    y = x .* (x>0) + 0.01 * x .* (x<=0);
end

function y = step(x)
    y = x > 0;
end