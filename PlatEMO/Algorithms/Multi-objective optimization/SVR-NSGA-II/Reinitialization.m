function Population=Reinitialization(Problem,NDS,ChangeCount)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Dec = NDS{ChangeCount+1}.decs;
    P   = 4;
    Up  = Problem.upper;
    Lp  = Problem.lower;
    if ChangeCount < P
        for i = 1 : size(Dec,1)
            for z = 1 : size(Dec,2)
                for j = 1 : ChangeCount+1
                    index     = NDS{j}.decs;
                    data(i,j) = index(i,z);
                end
                Dec(i,z) = modeltrain(data,Up(z),Lp(z));
            end
        end
    else
        for i = 1 : size(Dec,1)
            for z = 1 : size(Dec,2)
                for j = 1 : P+1
                    index     = NDS{ChangeCount-P+j}.decs;
                    data(i,j) = index(i,z);
                end
                Dec(i,z) = modeltrain(data,Up(z),Lp(z));
            end
        end
    end
    Population = Problem.Evaluation(Dec);
end

function predata = modeltrain(data,Up,Lp)
    x_train = data(:, 1:end-1);
    y_train = data(:, end);
    svr     = fitrsvm(x_train, y_train, 'KernelFunction', 'rbf', 'Epsilon', 0.05, 'BoxConstraint', 1e3);
    a       = predict(svr,data(:,2:end));
    predata = a(end);
    predata = max(min(predata, Up), Lp);
end