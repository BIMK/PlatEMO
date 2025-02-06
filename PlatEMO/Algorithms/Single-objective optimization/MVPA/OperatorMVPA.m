function Players_NEW = OperatorMVPA(Problem,Players,MVP,TEAMS,OBJS_TEAM,FranchisePlayer,npt,c1,c2)
%OperatorMVPA - The operator of Most Valuable Player Algorithm.
%
%   Players_NEW = OperatorMVPA(Problem,Players,MVP,TEAMS,OBJS_TEAM,FranchisePlayer,npt,c1,c2)
%
%   Example:
%       Players_NEW = OperatorMVPA(Problem,Players,MVP,TEAMS,OBJS_TEAM,FranchisePlayer,npt,c1,c2)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by H. R. E. H. Bouchekara (email: bouchekara.houssem@gmail.com)

    %% Parameter setting
    PlayersDec         = Players.decs;
    MVPDec             = MVP.decs;
    FranchisePlayerDec = FranchisePlayer.decs;
    [~,PROBLEMSIZE]    = size(PlayersDec);
    TeamsSize          = size(TEAMS,1);

    %% MVPA optimization
    for i = 1 : TeamsSize
        % Selecting the first team
        TEAMi = nonzeros(TEAMS(i,:));
        % Selecting the second team
        j = randi(TeamsSize); while j==i, j=randi(TeamsSize); end
        % Competition phase
        % Update of players skills of TEAMi
        PlayersDec(TEAMi,:) = PlayersDec(TEAMi,:)+...
            (c1*rand(npt(i),PROBLEMSIZE)).*(repmat(MVPDec,npt(i),1)-PlayersDec(TEAMi,:))+...
            (c2*rand(npt(i),PROBLEMSIZE)).*(repmat(FranchisePlayerDec(i,:),npt(i),1)-PlayersDec(TEAMi,:));

        TEAMWin = WinOrLoss(OBJS_TEAM,i,j); % check which team is winning

        if TEAMWin == 1 % TEAMi beats TEAMj
            PlayersDec(TEAMi,:) = PlayersDec(TEAMi,:)+rand*(PlayersDec(TEAMi,:)-repmat(FranchisePlayerDec(j,:),npt(i),1));
        else % TEAMj beats TEAMi
            PlayersDec(TEAMi,:) = PlayersDec(TEAMi,:)+rand*(repmat(FranchisePlayerDec(j,:),npt(i),1)-PlayersDec(TEAMi,:));
        end
        PlayersDec = bound(Problem,PlayersDec);

        Players_NEW = Problem.Evaluation(PlayersDec);
    end
end

function TEAMWin = WinOrLoss(Fitness,i,j)
    FitnessN = ((Fitness)-min(Fitness));
    p = 1-[FitnessN(i).^1/(sum(FitnessN(i).^1+FitnessN(j).^1)) FitnessN(j).^1/(sum(FitnessN(i).^1+FitnessN(j).^1))];
    if p(1) == p(2)
        if rand > 0.5
            TEAMWin = 1;
        else
            TEAMWin = 2;
        end
    else
        [p,ip]  = sort(p,'descend');
        TEAMWin = ip(1);
        if rand > p(1)
            TEAMWin = ip(2);
        end
    end
end

function x = bound(Problem,x)
    x = min(max(x,Problem.lower),Problem.upper);
end