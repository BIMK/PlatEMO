function DA = UpdateDA(DA,New,MaxSize,p,Problem,k)
% Update DA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Find the non-dominated solutions
    DA = [DA,New];
    ND = NDQSort(DA.objs,1);
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
        objmatrix = Popobj(~Choose,:);
        [~,Extreme(i)]     = min(max(Popobj(~Choose,:)./repmat(w(i,:),sum(~Choose),1),[],2)+0.1*objmatrix(:,i)/(1e-6)); % 带惩罚的AFS函数值
        Choose(Extreme(i)) = true; 
    end    
    
    %% Delete or add solutions to make a total of K solutions be chosen by truncation
    if sum(Choose) > MaxSize
        % Randomly delete several solutions
        Choosed = find(Choose);
        k       = randperm(sum(Choose),sum(Choose)-MaxSize);
        Choose(Choosed(k)) = false;
    elseif sum(Choose) < MaxSize
        % Calculate the angle
        angle    = acos(1-pdist2(Popobj,Popobj,'cosine'));
        Lp       = Shape_Estimate(DA,MaxSize,Zmin,Zmax);
        PopobjLP = DA.objs - Zmin;
        PopobjLP = PopobjLP+10^-6;
        tran_Obj = PopobjLP./repmat((sum(PopobjLP.^Lp,2)).^(1/Lp),1,M);
        for i = 1 : size(PopobjLP,1)
            DAPT(i) = norm(PopobjLP(i,:) - tran_Obj(i,:),p);
        end

        % Select the rest individuals
        while sum(Choose) < MaxSize
            if k == 0
                Select   = find(Choose);
                Remain   = find(~Choose);
                a        = min(angle(Remain,Select),[],2);
                bb       = a'.*(1 + Problem.FE^2*2/Problem.maxFE^2);
                cc       = 1./(DAPT(Remain).*(2-Problem.FE^2/Problem.maxFE^2));
                fitness  = bb.*(cc);
                [~,rho2] = max(fitness);
                Choose(Remain(rho2)) = true;
            elseif k == 1
                Select  = find(Choose);
                Remain  = find(~Choose);
                [~,rho] = max(min(angle(Remain,Select),[],2));
                Choose(Remain(rho)) = true;
            end
        end
    end
    DA = DA(Choose);
end