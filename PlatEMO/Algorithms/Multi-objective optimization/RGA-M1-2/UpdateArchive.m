function Archive = UpdateArchive(Population,N)
% Update Archive based on CDP

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
  
    %% Delete the duplicated solutions in current Population
    [~, Unduplicated] = unique(Population.objs,'rows');
    Population        = Population(Unduplicated);
    if length(Population) > N
        %% Non-dominated sorting    
        [FrontNo,MaxFNo] = NDSort(Population.objs,sum(max(0,Population.cons),2),N);
        Next = FrontNo < MaxFNo;
        
        %% Calculate the crowding distance of each solution
        CrowdDis = CrowdingDistance(Population.objs,FrontNo);
        
        %% Select the solutions in the last front based on their crowding distances
        Last     = find(FrontNo==MaxFNo);
        [~,Rank] = sort(CrowdDis(Last),'descend');
        Next(Last(Rank(1:N-sum(Next)))) = true;
        Archive  = Population(Next);
    else
        Archive = Population;
    end
end