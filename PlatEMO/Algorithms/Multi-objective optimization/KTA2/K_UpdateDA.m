function [DAobj,DAdec,DAvar] = K_UpdateDA(DAobj,DAdec,DAvar,MaxSize,p)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhenshou Song

    DA_Nor_pre = (DAobj - repmat(min(DAobj,[],1),size(DAobj,1),1))./repmat(max(DAobj,[],1) - min(DAobj,[],1),size(DAobj,1),1);

    %% Find the non-dominated solutions
    ND = NDSort(DAobj,1);
    DAobj = DAobj((ND==1),:);
    DAdec = DAdec((ND==1),:);
    DAvar = DAvar((ND==1),:);
    DA_Nor_pre = DA_Nor_pre((ND==1),:);
    N  = size(DAobj,1);
    if N <= MaxSize
        return;
    end
    
    %% Select the extreme solutions first
    Choose = false(1,N);
    
    M = size(DA_Nor_pre,2);
    select = randperm(M);
    Choose(select(1)) = true;

    %% Delete or add solutions to make a total of K solutions be chosen by truncation
    if sum(Choose) > MaxSize
        % Randomly delete several solutions
        Choosed = find(Choose);
        k = randperm(sum(Choose),sum(Choose)-MaxSize);
        Choose(Choosed(k)) = false;
    elseif sum(Choose) < MaxSize
        % Add several solutions by truncation strategy
        Distance = inf(N);
        for i = 1 : N-1
            for j = i+1 : N
                Distance(i,j) = norm(DA_Nor_pre(i,:)-DA_Nor_pre(j,:),p);
                Distance(j,i) = Distance(i,j);
            end
        end
        while sum(Choose) < MaxSize
            Remain = find(~Choose);
            [~,x]  = max(min(Distance(~Choose,Choose),[],2));
            Choose(Remain(x)) = true;
        end
    end
    DAobj = DAobj(Choose,:);
    DAdec = DAdec(Choose,:);
    DAvar = DAvar(Choose,:);
end