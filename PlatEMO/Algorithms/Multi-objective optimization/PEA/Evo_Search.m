function [PopDec,PopObj,PopCon,ObjMSE,ConMSE] = Evo_Search(P,wmax,Model_obj,Model_con,Problem)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhiyao Zhang (email: zhiyao_zhang0@163.com)

    PopDec = P.decs;
    PopObj = P.objs;
    PopCon = P.cons;
    ObjMSE = zeros(Problem.N,Problem.M);
    ConMSE = zeros(Problem.N,size(PopCon,2));
    w      = 1;
    while w <= wmax
        OffDec = OperatorGA(Problem,PopDec);
        [OffObj,Off_ObjMSE,OffCon,Off_ConMSE] = model_predict(Model_obj,Model_con,OffDec);

        PopDec = [PopDec;OffDec];
        PopObj = [PopObj;OffObj];
        PopCon = [PopCon;OffCon];
        ObjMSE = [ObjMSE;Off_ObjMSE];
        ConMSE = [ConMSE;Off_ConMSE];
        
        index  = EnvironmentalSelection(PopObj,ObjMSE,PopCon,ConMSE,length(P));
      
        PopDec = PopDec(index,:);
        PopObj = PopObj(index,:);
        PopCon = PopCon(index,:);
        ObjMSE = ObjMSE(index,:);
        ConMSE = ConMSE(index,:);
        w = w + 1;
    end
end

function [OffObj,Off_ObjMSE,OffCon,Off_ConMSE] = model_predict(Model_obj,Model_con,OffDec)
    N          = size(OffDec,1);
    Len_obj    = length(Model_obj);
    Len_con    = length(Model_con);
    OffObj     = zeros(N,Len_obj);
    OffCon     = zeros(N,Len_con);
    Off_ObjMSE = zeros(N,Len_obj);
    Off_ConMSE = zeros(N,Len_con);
      
    for i = 1 : N
        for j = 1 : Len_obj
            [OffObj(i,j),~,Off_ObjMSE(i,j)] = predictor(OffDec(i,:),Model_obj{j});
        end
        for j = 1 : Len_con
            [OffCon(i,j),~,Off_ConMSE(i,j)] = predictor(OffDec(i,:),Model_con{j});
        end
    end
end

function Next = EnvironmentalSelection(PopObj,ObjMSE,PopCon,ConMSE,N)
     %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort_CPPD(PopObj,ObjMSE,PopCon,ConMSE,N);
    Next = FrontNo < MaxFNo;
    Last = find(FrontNo == MaxFNo);

    %% Select the solutions in the last front
    if MaxFNo == 1
        Del  = Truncation(PopObj(Last,:),N);
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
    Next1 = 1:N1;
    Next2 = N1+1:N;
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
    [N,~]  = size(PopObj);
    
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