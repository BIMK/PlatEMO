function [Population,FrontNo,KneePoints] = EnvironmentalSelection(Population,FrontNo,MaxFNo,KneePoints,Distance,K)
% The environmental selection of KnEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Select the solutions in the first several fronts
    Next = FrontNo < MaxFNo;
    
    %% Select all the knee points in the last front
    Next(KneePoints) = true;
    
    %% Delete or add solutions to make a total of K solutions be chosen in the last front
    if sum(Next) < K
        Temp = find(FrontNo==MaxFNo & KneePoints==0);
        [~,Rank] = sort(Distance(Temp),'descend');
        Next(Temp(Rank(1:(K-sum(Next))))) = true;
    elseif sum(Next) > K
        Temp = find(FrontNo==MaxFNo & KneePoints==1);
        [~,Rank] = sort(Distance(Temp));
        Next(Temp(Rank(1:(sum(Next)-K)))) = false;
    end
    
    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    KneePoints = KneePoints(Next);
end