function [Population,FrontNo,App,Crowd] = EnvironmentalSelection(Population,P,theta,N,Center,R)
% The environmental selection of GFM-MOEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,N);
    Next             = find(FrontNo<=MaxFNo);

    %% Environmental selection
    [App,Dis] = subFitness(Population(Next).objs,P,Center,R);
    Choose    = LastSelection(Population(Next).objs,FrontNo(Next),App,Dis,theta,N);

    %% Population for next generation
    Population = Population(Next(Choose))  ;
    FrontNo    = FrontNo(Next(Choose));
    App        = App(Choose);
    Dis        = sort(Dis(Choose,Choose),2);
    Crowd      = Dis(:,1) + 0.1*Dis(:,2);
end

function Choose = LastSelection(PopObj,FrontNo,App,Dis,theta,N)
% Select part of the solutions in the last front

    %% Identify the extreme solutions
    NDS = find(FrontNo==1);
    [~,Extreme] = min(repmat(sqrt(sum(PopObj(NDS,:).^2,2)),1,size(PopObj,2)).*sqrt(1-(1-pdist2(PopObj(NDS,:),eye(size(PopObj,2)),'cosine')).^2),[],1);
    nonExtreme  = ~ismember(1:length(FrontNo),NDS(Extreme));
    %% Environmental selection
    Last   = FrontNo == max(FrontNo);
    Choose = true(1,size(PopObj,1));
    
    %% Non-dominated sort convergence and diversity
    while sum(Choose) > N
        Remain    = find(Choose&Last&nonExtreme);
        dis       = sort(Dis(Remain,Choose),2);
        dis       = dis(:,1) + 0.1*dis(:,2);
        fitness   = theta*dis + (1-theta)*App(Remain);
        [~,worst] = min(fitness);
        Choose(Remain(worst)) = false;
    end
end