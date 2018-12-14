function CA = UpdateCA(CA,New,MaxSize)
% Update CA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    CA = [CA,New];
    N  = length(CA);
    if N <= MaxSize
        return;
    end
    
    %% Calculate the fitness of each solution
    CAObj = CA.objs;
    CAObj = (CAObj-repmat(min(CAObj),N,1))./(repmat(max(CAObj)-min(CAObj),N,1));
    I = zeros(N);
    for i = 1 : N
        for j = 1 : N
            I(i,j) = max(CAObj(i,:)-CAObj(j,:));
        end
    end
    C = max(abs(I));
    F = sum(-exp(-I./repmat(C,N,1)/0.05)) + 1;
    
    %% Delete part of the solutions by their fitnesses
    Choose = 1 : N;
    while length(Choose) > MaxSize
        [~,x] = min(F(Choose));
        F = F + exp(-I(Choose(x),:)/C(Choose(x))/0.05);
        Choose(x) = [];
    end
    CA = CA(Choose);
end