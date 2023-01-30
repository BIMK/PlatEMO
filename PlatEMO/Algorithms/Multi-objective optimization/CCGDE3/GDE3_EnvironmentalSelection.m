function Population = GDE3_EnvironmentalSelection(Population,Offspring,N)
% The environmental selection of CCGDE3

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    % Select the offsprings dominating their corresponding parents firstly,
    % and then select the offspring population non-dominated with the
    % parent population
    PopObj = Population.objs;
    OffObj = Offspring.objs;
    % The offsprings which can replace its parent
    updated = all(PopObj>=OffObj,2);
    % The offsprings which can add to the population
    selected = any(PopObj<OffObj,2) & any(PopObj>OffObj,2);
    % Update the population
    Population(updated) = Offspring(updated);
    Population          = [Population,Offspring(selected)];
    
    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,N);
    Next = FrontNo < MaxFNo;
    
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(Population.objs,FrontNo);
    
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% Population for next generation
    Population = Population(Next);
end