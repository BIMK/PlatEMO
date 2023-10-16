function [Population,FrontNo,CrowdDis] = EnvironmentalSelection(Problem,Population,N,type,eta)
% Environmental selection

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Non-dominated sorting
    PopObj = Population.objs;
    if type
        % Selection of DV1
        for i = 1 : length(Population)
            PopX         = Problem.Perturb(Population(i).dec);
            PopObjV(i,:) = mean(PopX.objs,1);
        end
        DR = sum(abs(PopObjV-PopObj),2);
        DR(DR<=eta) = 0;
        [FrontNo,MaxFNo] = NDSort(PopObj,DR,N);
    else
        % Selection of DV2
        [FrontNo,MaxFNo] = NDSort(PopObj,N);
    end
    Next = FrontNo < MaxFNo;
    
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(PopObj,FrontNo);
    
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next); 
end