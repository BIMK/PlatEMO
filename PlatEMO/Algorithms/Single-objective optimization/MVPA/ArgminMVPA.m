function [MVP,TEAMS,OBJS_TEAM,FranchisePlayer]=ArgminMVPA(Players,TEAMS)
%ArgminMVPA - Update MVp and FranchisePlayer for the Most Valuable Player Algorithm.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by H. R. E. H. Bouchekara (email: bouchekara.houssem@gmail.com)

    TeamsSize = size(TEAMS,1);
    OBJS      = Players.objs;
    for i = 1 : TeamsSize
        [~,~,nTEAM] = find(TEAMS(i,:));
        [OBJS_TEAM(i),j] = min(OBJS(nTEAM));
        FranchisePlayer_ID(i) = nTEAM(j);
    end
    [~,MVP_ID] = min(OBJS);
    FranchisePlayer = Players(FranchisePlayer_ID);
    MVP = Players(MVP_ID);
end