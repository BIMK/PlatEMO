function Next = Pop_Update(PopObj,PopCon,N)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhiyao Zhang (email: zhiyao_zhang0@163.com)

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(PopObj,PopCon,N);
    Next = FrontNo < MaxFNo;
    Last = find(FrontNo == MaxFNo);
    
    %% Select the solutions in the last front
    zmin   = min(PopObj,[],1); zmax = max(PopObj,[],1);
    PopObj = (PopObj - zmin)./max(zmax - zmin,10e-10);
    PopCon = max(PopCon,0);
    PopCon = PopCon./max(PopCon,[],1);
    CV     = sum(PopCon,2);
    CV     = CV./max(max(CV,10e-10));
    
    if MaxFNo == 1
        Del  = Truncation(PopObj(Last,:),CV(Last,:),N);
        Next(Last(Del)) = true; 
    else
        Choose = Dist_Selection(PopObj(Next,:),CV(Next,:),PopObj(Last,:),CV(Last,:),N - sum(Next));
        Next(Last(Choose)) = true;
    end
end

function Choose = Dist_Selection(PopObj1,CV1,PopObj2,CV2,mu)
    PopObj   = [PopObj1,CV1;PopObj2,CV2];
    N        = size(PopObj,1);
    N1       = size(PopObj1,1);
    Distance = inf(N);
    
    %% Calculate the shifted distance between each two solutions
    for i = 1 : N
        SPopObj = max(PopObj,repmat(PopObj(i,:),N,1));
        for j = 1 : N
            Distance(i,j) = norm(PopObj(i,:)-SPopObj(j,:),2);
        end
    end
    Distance(logical(eye(length(Distance)))) = inf;
    
    %% Calculate D
    Next1 = 1 : N1;
    Next2 = N1+1 : N;
    for i = 1 : mu
        Distance1 = sort(Distance(Next2,Next1),2);
        [~,index] = max(Distance1(:,1));
        Next1     = [Next1,Next2(index)];
        Next2(index) = [];
    end
    Choose = Next1(N1+1:end) - N1;
end

function Del = Truncation(PopObj,CV,K)
    %% Select part of the solutions by truncation
    [N,~]  = size(PopObj);
    PopObj = [PopObj,CV];
    
    %% Calculate the shifted distance between each two solutions
    Distance = inf(N);
    for i = 1 : N
        SPopObj = max(PopObj,repmat(PopObj(i,:),N,1));
        for j = 1 : N
            Distance(i,j) = norm(PopObj(i,:)-SPopObj(j,:),2);
        end
    end
    
    %% Truncation
    Distance(logical(eye(length(Distance)))) = inf;
    Del = true(1,N);
    while sum(Del) > K
        Remain   = find(Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = false;
    end
end