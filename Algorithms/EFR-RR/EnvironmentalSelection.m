function [Population,z,znad] = EnvironmentalSelection(Population,W,N,K,z,znad)
% The environmental selection of EFR-RR

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Normalization
    [PopObj,z,znad] = Normalization(Population.objs,z,znad);

    %% Environmental selection
    RgFrontNo = MaximumRanking(PopObj,W,K);
    MaxFNo    = find(cumsum(hist(RgFrontNo,1:max(RgFrontNo)))>=N,1);
    LastFront = find(RgFrontNo==MaxFNo);
    LastFront = LastFront(randperm(length(LastFront)));
    RgFrontNo(LastFront(1:sum(RgFrontNo<=MaxFNo)-N)) = inf;
    Next      = RgFrontNo <= MaxFNo;
    % Population for next generation
    Population = Population(Next);
end