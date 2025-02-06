function [Population,z,znad,z_c,znad_c] = EnvironmentalSelection(Population,W,N,z,znad,z_c,znad_c)
% The environmental selection of theta-DEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming (email: 20151000334@cug.edu.cn)

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,N);
    St = find(FrontNo<=MaxFNo);

    %% Normalization
    [PopObj,PopCon,z,znad,z_c,znad_c] = Normalization(Population(St).objs,Population(St).cons,z,znad,z_c,znad_c);
    CV      = sum(max(0,PopCon),2);
    fr      = sum(CV==0)/N;
  
    %% theta-non-dominated sorting
    tFrontNo = tNCDSort(PopObj,PopCon,W,fr);
    MaxFNo    = find(cumsum(hist(tFrontNo,1:max(tFrontNo)))>=N,1);
    LastFront = find(tFrontNo==MaxFNo);
    LastFront = LastFront(randperm(length(LastFront)));
    tFrontNo(LastFront(1:sum(tFrontNo<=MaxFNo)-N)) = inf;
    Next      = St(tFrontNo<=MaxFNo);
    % Population for next generation
    Population = Population(Next);
end