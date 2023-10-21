function DA = UpdateDA(DA,New,MaxSize)
% Update DA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Find the non-dominated solutions
    DA = [DA,New];
    ND = NDSort(DA.objs,1);
    DA = DA(ND==1);
    N  = length(DA);
    if N <= MaxSize
        return;
    end
    Popobj = DA.objs;
    
    %% Normalization
    Zmin   = min(Popobj,[],1);
    Zmax   = max(Popobj,[],1);
    Popobj = (Popobj-repmat(Zmin,size(Popobj,1),1))./repmat(Zmax-Zmin,size(Popobj,1),1);
    Choose = false(1,N);

    %% Slect the extreme solutions first
    M = size(Popobj,2);    
    Extreme = zeros(1,M);
    w       = zeros(M)+1e-6+eye(M);
    for i = 1 : M
        objmatrix = DA(~Choose).objs;
        [~,Extreme(i)] = min(max(DA(~Choose).objs./repmat(w(i,:),sum(~Choose),1),[],2)+0.1*objmatrix(:,i)/(1e-6)); % 带惩罚的AFS函数值
        Choose(Extreme(i)) = true; 
    end    

    %% Delete or add solutions to make a total of K solutions be chosen by truncation
    if sum(Choose) > MaxSize
        % Randomly delete several solutions
        Choosed = find(Choose);
        k = randperm(sum(Choose),sum(Choose)-MaxSize);
        Choose(Choosed(k)) = false;
    elseif sum(Choose) < MaxSize
        % Calculate the angle
        angle = acos(1-pdist2(Popobj,Popobj,'cosine'));
        % Select the rest individuals
        while sum(Choose) < MaxSize
            Select  = find(Choose);
            Remain  = find(~Choose);
            [~,rho] = max(min(angle(Remain,Select),[],2));
            Choose(Remain(rho)) = true;
        end
    end
    DA = DA(Choose);
end