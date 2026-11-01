function NewCluster = ClusterEnvironmentalSelection(Cluster,Offspring1,W,Z)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yubo Wang (email: 1202111230@cug.edu.cn)

    numClusters     = size(W,1);
    NewCluster      = cell(numClusters,1);
    TSwarm          = [Cluster{:}];
    TSwarm          = [TSwarm Offspring1];
    [N,~]           = size(TSwarm.objs);
    Fitness         = CalFitness(TSwarm.objs,TSwarm.cons,0);
    [~,Association] = min(acos(1-pdist2((TSwarm.objs-repmat(Z,N,1)),W,'cosine')),[],2);
    for i = 1 : numClusters
        ParentNum = length(Cluster{i});
        FinalNum  = ParentNum/2;
        Flag      = (Association==i);
        InSwarm   = TSwarm(Flag);
        InFit     = Fitness(Flag);
        OutSwarm  = TSwarm(~Flag);
        if length(InSwarm) >= FinalNum
            [~,Rank]      = sort(InFit);
            NewCluster{i} = [InSwarm(Rank(1:FinalNum))];
        else
            Theta = acos(1-pdist2((OutSwarm.objs-repmat(Z,length(OutSwarm),1)),W(i,:),'cosine'));
            [~,Rank]      = sort(Theta);
            NewCluster{i} = [InSwarm OutSwarm(Rank(1:(FinalNum-length(InSwarm))))];
        end
    end
end