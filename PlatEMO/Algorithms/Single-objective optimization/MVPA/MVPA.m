classdef MVPA < ALGORITHM
% <2020> <single> <real/integer> <large/none> <constrained/none>
% Most valuable player algorithm

%------------------------------- Reference --------------------------------
% H. Bouchekara. Most valuable player algorithm: a novel optimization
% algorithm inspired from sport. Operational Research, 2020, 20: 139-195.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by H. R. E. H. Bouchekara (email: bouchekara.houssem@gmail.com)

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [TeamsSize,ArchiveSize,c1,c2] = Algorithm.ParameterSet(round(Problem.N/5),round(Problem.N/3),1,2);
            
            %% Generate random population
            Players = Problem.Initialization();
            [MVP,TEAMS,OBJS_TEAM,FranchisePlayer,npt] = InitializeMVPA(Problem,Players,TeamsSize);

            %% Optimization
            while Algorithm.NotTerminated(Players)
                [~,rank]         = sort(Players.objs);
                Archive          = Players(rank(1:ArchiveSize));
                Players_NEW      = OperatorMVPA(Problem,Players,MVP,TEAMS,OBJS_TEAM,FranchisePlayer,npt,c1,c2);
                replace          = find(Players_NEW.objs<=Players.objs);
                Players(replace) = Players_NEW(replace);
                [MVP,TEAMS,OBJS_TEAM,FranchisePlayer]=ArgminMVPA(Players,TEAMS);
                Players_All      = [Archive, Players];
                [~,rank1]        = sort(Players_All);
                Players          = Players_All(rank1);
            end
        end
    end
end