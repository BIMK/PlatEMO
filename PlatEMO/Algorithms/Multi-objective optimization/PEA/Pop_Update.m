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
    if MaxFNo == 1
        Del = Truncation(PopObj(Last,:),N);
        Next(Last(Del)) = true; 
    else
        Choose = Dist_Selection(PopObj(Next,:),PopObj(Last,:),N - sum(Next));
        Next(Last(Choose)) = true;
    end
end

function Choose = Dist_Selection(PopObj1,PopObj2,mu)
    PopObj = [PopObj1;PopObj2];
    N      = size(PopObj,1);
    N1     = size(PopObj1,1);
    
    %% Calculate the distance between each two solutions
    for i = 1 : N
        for j = i+1 : N
            Distance(i,j) = norm(PopObj(i,:)-PopObj(j,:),2);
            Distance(j,i) = Distance(i,j);
        end
    end
    
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

function Del = Truncation(PopObj,K)
    %% Select part of the solutions by truncation
    [N,~] = size(PopObj);

    %% Calculate the distance between each two solutions
    Distance = inf(N);
    for i = 1 : N
         for j = 1 : N
            Distance(i,j) = norm(PopObj(i,:) - PopObj(j,:),2);
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