function [Population,Dec,Mask,FrontNo,CrowdDis,sRatio] = EnvironmentalSelection(Population,Dec,Mask,N,len,num)
% The environmental selection of MOEA/PSL

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Delete duplicated solutions
    success = false(1,length(Population));
    [~,uni] = unique(Population.objs,'rows');
    if length(uni) == 1
        [~,uni] = unique(Population.decs,'rows');
    end
    Population = Population(uni);
    Dec        = Dec(uni,:);
    Mask       = Mask(uni,:);
    N          = min(N,length(Population));

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,N);
    Next = FrontNo < MaxFNo;    
    
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(Population.objs,FrontNo);
    
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% Calculate the ratio of successful offsprings
    success(uni(Next)) = true;
    s1     = sum(success(len+1:len+num));
    s2     = sum(success(len+num+1:end));
    sRatio = (s1+1e-6)./(s1+s2+1e-6);
    sRatio = min(max(sRatio,0.1),0.9);
    
    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
    Dec        = Dec(Next,:);
    Mask       = Mask(Next,:);
end