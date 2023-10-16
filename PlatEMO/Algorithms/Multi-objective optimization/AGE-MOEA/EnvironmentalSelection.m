function [Population,FrontNo,CrowdDis] = EnvironmentalSelection(Population,N)
% The environmental selection of AGE-MOEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Annibale Panichella

    %% let's round the objective values
    objs = round(Population.objs, 6);

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(objs,Population.cons,N);
    Next = FrontNo < MaxFNo;

    [nInd, ~]    = size(objs);
    CrowdDis = zeros(1,nInd);
    
    %% Calculate the crowding distance of each solution
    front1 = objs(FrontNo==1,:);
    if (size(front1,1)>1)
        IdealPoint = min(front1);
    else
        IdealPoint = front1;
    end
        
    [CrowdDis(FrontNo==1), p, normalization] = SurvivalScore(front1, IdealPoint);  
    for i=2:MaxFNo
        front = objs(FrontNo==i,:);
        [m,~] = size(front);
        front = front./repmat(normalization',m,1);
        CrowdDis(FrontNo==i) = 1./pdist2(front, IdealPoint, 'minkowski', p);  
    end
    
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
end