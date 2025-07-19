function Dec = MLP(DecSource,ChangeCount,P,score)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Dec = DecSource{ChangeCount+1};
    
    if ChangeCount<P
        for i = 1 : size(Dec,1)
            for j = 1 : ChangeCount+1
                index     = DecSource{j};
                data(j,:) = index(i,:);
            end
            Dec(i,:) = modeltrain(data,score);
        end
    else
        for i = 1 : size(Dec,1)
            for j = 1 : P+1
                index     = DecSource{ChangeCount-P+j};
                data(j,:) = index(i,:);
            end
            Dec(i,:) = modeltrain(data,score);
        end
    end
end

function predata = modeltrain(data,score)
    PO      = find(score~=0);
    predata = zeros(1,size(data,2));
    
    data        = data(:,PO);
    train_data  = data(1:end-1,:);
    target_data = data(end,:);

    input_size  = size(train_data, 2);
    hidden_size = 10;
    output_size = size(target_data, 2);
    W1 = randn(input_size, hidden_size);
    W2 = randn(hidden_size, output_size);
    b1 = randn(1, hidden_size);
    b2 = randn(1, output_size);

    learning_rate = 0.1;
    num_epochs    = 500;

    for epoch = 1 : num_epochs
        hidden_output = sigmoid(train_data* W1 + b1);
        output        = sigmoid(hidden_output * W2 + b2);

        d_output = output - target_data;
        d_hidden_output = d_output * W2' .* hidden_output .* (1 - hidden_output);
        d_W2 = hidden_output' * d_output;
        d_b2 = sum(d_output);
        d_W1 = train_data' * d_hidden_output;
        d_b1 = sum(d_hidden_output);

        W1 = W1 - learning_rate * d_W1;
        W2 = W2 - learning_rate * d_W2;
        b1 = b1 - learning_rate * d_b1;
        b2 = b2 - learning_rate * d_b2;
    end
    input_data    = target_data ;
    hidden_output = sigmoid(input_data * W1 + b1);
    predata(:,PO) = sigmoid(hidden_output * W2 + b2);
end

function y = sigmoid(x)
    y = 1./(1 + exp(-x));
end