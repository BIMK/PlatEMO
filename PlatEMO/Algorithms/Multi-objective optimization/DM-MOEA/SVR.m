function Mask = SVR(MaskSource,ChangeCount,P)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Mask = MaskSource{ChangeCount+1};
    data = zeros(P+1,size(Mask,2));
    for i = 1 : size(Mask,1)
        for j = 1 : P+1
            index     = MaskSource{ChangeCount-P+j};
            data(j,:) = index(i,:);
        end
        Mask(i,:) = train(data);
    end
end

function predata = train(data)
    predata = zeros(1, size(data,2));
    score   = mean(data,1);
    index1  = find(score==1);
    index2  = find(score==0);
    index3  = setdiff(1:size(data,2),[index1,index2]);
    data    = data(:,index3);
    [N,D]   = size(data);
    X_train = data(1:N-1, :);
    Y_train = data(2:N, :);
    X_test  = data(N, :);

    X_train = double(X_train);
    Y_train = double(Y_train);
    X_test  = double(X_test);

    y_pred = zeros(1, D);

    for i = 1 : D
        svm = fitrsvm(X_train, Y_train(:, i), 'KernelFunction','RBF');
        y_pred(i) = predict(svm, X_test)>rand;
    end
    if ~isempty(index1)
        predata(:,index1) = 1;
    end
    predata(:,index2) = 0;
    predata(:,index3) = y_pred;
end