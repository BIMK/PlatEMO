function [Population,FrontNo,SpCrowdDis] = non_domination_scd_sort(Population,N)
% Sort the population according to non-dominated relationship and special crowding 

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    %% Non-dominated sorting
    N = min(N,length(Population));
    [FrontNo,MaxFNo] = NDSort(Population.objs,N);
    Next = FrontNo < MaxFNo;
    
    %% Calculate the special crowding distance of each solution
    SpCrowdDis_Obj = ModifiedCrowdingDistance(Population.objs,FrontNo);
    SpCrowdDis_Dec = ModifiedCrowdingDistance(Population.decs,FrontNo);
    SpCrowdDis     = max(SpCrowdDis_Obj,SpCrowdDis_Dec);
    for i = 1 : MaxFNo
        Front   = find(FrontNo==i);
        Avg_Obj = mean(SpCrowdDis_Obj(Front));
        Avg_Dec = mean(SpCrowdDis_Dec(Front));
        replace = SpCrowdDis_Obj(Front)<=Avg_Obj & SpCrowdDis_Dec(Front)<=Avg_Dec;
        SpCrowdDis(Front(replace)) = min(SpCrowdDis_Obj(Front(replace)),SpCrowdDis_Dec(Front(replace)));
    end
    
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(SpCrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% Limit the size of Population 
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    SpCrowdDis = SpCrowdDis(Next);
end

function CrowdDis = ModifiedCrowdingDistance(PopObj,FrontNo)
    [N,M]    = size(PopObj);
    CrowdDis = zeros(1,N);
    Fronts   = setdiff(unique(FrontNo),inf);
    for f = 1 : length(Fronts)
        Front = find(FrontNo==Fronts(f));
        Fmax  = max(PopObj(Front,:),[],1);
        Fmin  = min(PopObj(Front,:),[],1);
        for i = 1 : M
            [~,Rank] = sortrows(PopObj(Front,i));
            CrowdDis(Front(Rank(1))) = CrowdDis(Front(Rank(1))) + 1;
            for j = 2 : length(Front)-1
                if Fmax(i) == Fmin(i)
                    CrowdDis(Front(Rank(j))) = CrowdDis(Front(Rank(j)))+1;
                else
                    CrowdDis(Front(Rank(j))) = CrowdDis(Front(Rank(j)))+(PopObj(Front(Rank(j+1)),i)-PopObj(Front(Rank(j-1)),i))/(Fmax(i)-Fmin(i));
                end
            end
        end
    end
end