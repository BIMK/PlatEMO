function Population = EnvironmentalSelection(Population,Offspring,N)
% The environmental selection of GDE3

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    %% Select by constraint-domination
    PopObj    = Population.objs;
    PopCon    = Population.cons;
    feasibleP = all(PopCon<=0,2);
    OffObj    = Offspring.objs;
    OffCon    = Offspring.cons;
    feasibleO = all(OffCon<=0,2);
    % The offsprings which can replace its parent
    updated = ~feasibleP&feasibleO  | ...
              ~feasibleP&~feasibleO & all(PopCon>=OffCon,2) | ...
              feasibleP&feasibleO   & all(PopObj>=OffObj,2);
    % The offsprings which can add to the population
    selected = feasibleP&feasibleO & any(PopObj<OffObj,2) & any(PopObj>OffObj,2);
    % Update the population
    Population(updated) = Offspring(updated);
    Population          = [Population,Offspring(selected)];
    
    %% Select by non-dominated sorting and crowding distance
    PopObj   = Population.objs;
    PopCon   = Population.cons;
    feasible = all(PopCon<=0,2);
    % Non-dominated sorting based on constraint-domination
    FrontNo = inf(1,length(Population));
    [FrontNo(feasible),MaxFNo] = NDSort(PopObj(feasible,:),inf);
    FrontNo(~feasible) = NDSort(PopCon(~feasible,:),inf) + MaxFNo;
    % Determine the last front
    MaxFNo    = find(cumsum(hist(FrontNo,1:max(FrontNo)))>=N,1);
    lastFront = find(FrontNo==MaxFNo);
    % Eliminate solutions in the last front one by one
    while length(lastFront) > N - sum(FrontNo<MaxFNo)
        [~,worst] = min(CrowdingDistance(PopObj(lastFront,:)));
        lastFront(worst) = [];
    end
    Population = Population([find(FrontNo<MaxFNo),lastFront]);
end

function CrowdDis = CrowdingDistance(PopObj)
% Calculate the crowding distance of each solution in the same front

    [N,M]    = size(PopObj);
    
    CrowdDis = zeros(1,N);
    Fmax     = max(PopObj,[],1);
    Fmin     = min(PopObj,[],1);
    for i = 1 : M
        [~,rank] = sortrows(PopObj(:,i));
        CrowdDis(rank(1))   = inf;
        CrowdDis(rank(end)) = inf;
        for j = 2 : N-1
            CrowdDis(rank(j)) = CrowdDis(rank(j))+(PopObj(rank(j+1),i)-PopObj(rank(j-1),i))/(Fmax(i)-Fmin(i));
        end
    end
end