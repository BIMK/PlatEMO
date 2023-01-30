function [AA,gBest] = A_get(Problem,Population,A,iter)
% The updating strategy of archive(A) in S-ECSO

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    x = Population.decs;
    xObj_fitness = Population.objs;
    if iter == 1
        A       = [];
        Adecs   = [];
        value_A = [];
    else
        Adecs   = A.decs;
        value_A = A.objs;
    end

    %% update A
    Compare_Obj = [value_A;xObj_fitness];
    compare_x   = [Adecs;x];
    CC          = [A,Population];

    [FrontNo,~]   = NDSort(Compare_Obj,Inf);
    FrontNo_index = find(FrontNo == 1);

    Adecs   = compare_x(FrontNo_index,:);
    value_A = Compare_Obj(FrontNo_index,:);
    AA      = CC(FrontNo_index);


    [Adecs,index] = unique(Adecs,'rows');
    AA            = AA(index);

    %% adding perturbation
    ELS_A = AA.decs;
    for i = 1 : size(Adecs,1)
        j = ceil(rand*size(Adecs,2));
        ELS_A(i,j) = ELS_A(i,j)+ (Problem.upper(j)-Problem.lower(j))*normrnd(0,1);
        ELS_A(i,j) = min(max(ELS_A(i,j),Problem.lower(j)),Problem.upper(j));
    end

    ELS_A1 = Problem.Evaluation(ELS_A);

    %%  truncating (using SPEA2)
    AA    = [AA,ELS_A1];
    row_A = size(AA.decs,1);
    if row_A > Problem.N
        [FrontNo,MaxFNo] = NDSort(AA.objs,Problem.N);
        Next = false(1,length(FrontNo));
        Next(FrontNo<MaxFNo) = true;
        PopObj = AA.objs;
        fmax   = max(PopObj(FrontNo==1,:),[],1);
        fmin   = min(PopObj(FrontNo==1,:),[],1);
        PopObj = (PopObj-repmat(fmin,size(PopObj,1),1))./repmat(fmax-fmin,size(PopObj,1),1);

        %% Select the solutions in the last front
        Last = find(FrontNo==MaxFNo);
        del  = Truncation(PopObj(Last,:),length(Last)-Problem.N+sum(Next));
        Next(Last(~del)) = true;
        AA = AA(Next);
    end

    %% gBest
    [M,~]     = size(AA.decs);
    r         = rand(1,size(AA.objs,2));
    r_matr    = repmat(r,M,1);
    f_gBest   = sum(r_matr.*AA.objs,2)/sum(r);
    [~,index] = min(f_gBest);
    B         = AA.decs;
    gBest     = B(index,:);
end

function Del = Truncation(PopObj,K)
% Select part of the solutions by truncation
 
    N = size(PopObj,1);

    %% Truncation
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,N);
    while sum(Del) < K 
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end