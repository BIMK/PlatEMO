function Next = Pop_Reselect(PopObj,PopCon,N)

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
    zmin   = min(PopObj);
    zmax   = max(PopObj);
    PopObj = (PopObj - zmin)./max(zmax - zmin,10e-10);
    global phase
    if phase == 2
        [FrontNo,MaxFNo] = NDSort(PopObj,PopCon,N);
    else
        [FrontNo,MaxFNo] = NDSort(PopObj,N);
    end
    Next = FrontNo < MaxFNo;
    Last = find(FrontNo == MaxFNo);

    %% Select the solutions in the last front
    if MaxFNo == 1
        Del = Truncation(PopObj(Last,:),N);
        Next(Last(Del)) = true;
    else
        Choose = Dis_Selection(PopObj,Last,N-sum(Next));  
        Next(Last(Choose)) = true;
    end
end

function Choose = Dis_Selection(PopObj,Last,mu)
    N = size(PopObj,1);
    Distance = inf(N);
    for i = 1 : N
        for j = [1:i-1,i+1:N]
           Distance(i,j) = norm(PopObj(i,:)-PopObj(j,:),2);
        end
    end
    Distance  = sort(Distance,2);
    D         = 1./(Distance(:,1)+2);
    D         = D(Last);
    [~,index] = sort(D);
    Choose    = index(1:mu);
end

function Del = Truncation(PopObj,K)
    %% Select part of the solutions by truncation
    N = size(PopObj,1);
    
    %% Calculate the distance between each two solutions
    Distance = inf(N);
    for i = 1 : N
        for j = 1 : N
            Distance(i,j) = norm(PopObj(i,:)-PopObj(j,:),2);
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