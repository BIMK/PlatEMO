function [Population,Fitness,Next] = EnvironmentalSelectionT2(Population,N,alpha,gamma,para)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Mengjun Ming

    %% Parameter
    popSize   = length(Population);
    NumSeq    = (1:popSize)';
    RankConvg = zeros(popSize,1);
    RankDivs  = zeros(popSize,1);
    
    %% Modify the infeasible solutions
    PopObj = Population.objs;
    PopCon = Population.cons;
    z      = min(PopObj,[],1);
    n      = max(PopObj,[],1);
    
    Infeasible_all = any(PopCon>0,2);
    phi_max = max(sum(max(0,PopCon(Infeasible_all,:)),2));
    
    M          = length(z);
    [W,~]      = UniformPoint(N,M);
    [~,Region] = min(pdist2(PopObj-z,W,'cosine'),[],2);  
    PopObj_2   = PopObj;
    for i = 1:size(W,1)
        index = find(Region==i);
        if (~isempty(index))
            Objs_temp  = PopObj_2(index,:);
            Cons_temp  = PopCon(index,:);
            Infeasible = any(Cons_temp>0,2);
            if (sum(Infeasible)~=0)
                F_max = max(Objs_temp,[],1);
                PopObj_2(index(Infeasible),:) = Objs_temp(Infeasible,:)+(sum(max(0,Cons_temp(Infeasible,:)),2)/phi_max).^(exp(para)/max(gamma,0.000001)).*(F_max-Objs_temp(Infeasible,:));
            end
        end
    end

    %% Non-dominated sorting
    Dominate = false(popSize);
    for i = 1 : popSize-1
        for j = i+1 : popSize
            k = any(PopObj_2(i,:)<PopObj_2(j,:)) - any(PopObj_2(i,:)>PopObj_2(j,:));
            if k == 1
                Dominate(i,j) = true;
            elseif k == -1
                Dominate(j,i) = true;
            end
        end
    end
    
    %% Calculate S(i)
    S = sum(Dominate,2);
    
    %% Calculate R(i)
    R = zeros(1,popSize);
    for i = 1 : popSize
        R(i) = sum(S(Dominate(:,i)));
    end
    FrontNo = R + 1;
    
    %% Calculate the crowding distance of each solution
    PopObj   = Population.objs;
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Distance = sort(Distance,2);
    CrowdDis = Distance(:,floor(sqrt(popSize)));

    %% Add a middle column
    MiddleLevel = zeros(popSize,1);
    
    %% Environmental selection
    Next1 = FrontNo == 1;
    if sum(Next1) <= N
        [~,indx_Convg] = sortrows([FrontNo',-CrowdDis]);        
    elseif sum(Next1) > N
        Del  = Truncation(Population(Next1).objs,sum(Next1)-N);
        Temp = find(Next1);
        Next1(Temp(Del)) = false;
        MiddleLevel(Temp(Del)) = 1;
        [~,indx_Convg] = sortrows([FrontNo',MiddleLevel,-CrowdDis]);
    end
    RankConvg(indx_Convg) = NumSeq;
    
    %% Environmental selection -- diversity
    FrontNo_D = ones(popSize,1);
    for i = 1:size(W,1)
        index = find(Region==i);
        if (~isempty(index))
            Objs_temp = PopObj_2(index,:);          
            g_temp = sum((Objs_temp-z).*W(i,:),2);
            [~,index_FrontNo_D] = sort(g_temp);
            FrontNo_D(index(index_FrontNo_D)) = (1:length(g_temp))';
        end
    end
    [~,indx_divs] = sortrows([FrontNo_D,-CrowdDis]);
    RankDivs(indx_divs) = NumSeq;
    
    %% Population for next generation
    RankSolution = alpha*RankConvg+(1-alpha)*RankDivs;
    [~,Rank]     = sort(RankSolution);
    Next = zeros(1,length(Population));
    Population   = Population(Rank(1:N));
    Fitness = RankSolution(Rank(1:N));
    Choosed = Rank(1:N);
    Next(Choosed) = 1;
end

function Del = Truncation(PopObj,K)
% Select part of the solutions by truncation

    %% Truncation
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,size(PopObj,1));
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end