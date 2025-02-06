function [Population2,Fitness,Population3,CV] = EvenSearch(Population,non_dom,N,tau,max_cv)
% Even search for helper population

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming (email: 20151000334@cug.edu.cn)

    %% Determine if in promising region (1 is in)
    lables = GetLable(Population.objs,non_dom.objs);
    Next1 = lables == 1;
    First_Population = Population(Next1);
    
    % first choose non-dominated solutions near UPF
    Fitness = CalFitness(First_Population.objs);
    Next2 = Fitness < 1;
    if sum(Next2) > N
        Del  = Truncation(First_Population(Next2).objs,sum(Next2)-N);
        Temp = find(Next2);
        Next2(Temp(Del)) = false;
    end
    % store non-dominated solutions in former N
    Population2 = First_Population(Next2);
    Fitness = Fitness(Next2);
    
    % then choose CV better solutions near or in feasible regions
    % first limit the position in the objective space
    Second_Population = First_Population;
    CV = sum(max(0,Second_Population.cons),2);
    epsilon = (max_cv)*(1-tau)^2;
    Next3 = CV <= epsilon;
    if sum(Next3) > N
        Del  = Truncation(Second_Population(Next3).objs,sum(Next3)-N);
        Temp = find(Next3);
        Next3(Temp(Del)) = false;
    end
    % store more feasible solutions in former N
    Population3 = Second_Population(Next3);
    CV = CV(Next3);
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