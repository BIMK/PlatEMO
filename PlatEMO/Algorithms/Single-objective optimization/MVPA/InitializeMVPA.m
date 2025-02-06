function [MVP,TEAMS,OBJS_TEAM,FranchisePlayer,npt]=InitializeMVPA(Problem,Players,TeamsSize)
% The extra initialization for the Most Valuable Player Algorithm

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by H. R. E. H. Bouchekara (email: bouchekara.houssem@gmail.com)

    PlayersSize = Problem.N;
    [TEAMS,npt] = TeamsFormation(PlayersSize,TeamsSize);
    OBJS        = Players.objs;

    for i = 1 : TeamsSize
        [~,~,nTEAM] = find(TEAMS(i,:));
        [OBJS_TEAM(i),j] = min(OBJS(nTEAM));
        FranchisePlayer_ID(i) = nTEAM(j);
    end
    [~,MVP_ID] = min(OBJS);
    FranchisePlayer = Players(FranchisePlayer_ID);
    MVP = Players(MVP_ID);
end

function [TEAMS,npt] = TeamsFormation(PlayersSize,TeamsSize)
    Players_ID = 1 : PlayersSize;
    np1 = ceil(PlayersSize/TeamsSize);
    np2 = np1*TeamsSize-PlayersSize;
    % number of players per team
    npt = [repmat(np1,1,TeamsSize-np2) repmat(np1-1,1,np2)];
    for i = 1 : TeamsSize
        for j = 1 : npt(i)
            a = randi(length(Players_ID));
            TEAMS(i,j) = Players_ID(a);
            Players_ID(a) = [];
        end
        TEAMS(i,:) = sort(TEAMS(i,:));
    end
end