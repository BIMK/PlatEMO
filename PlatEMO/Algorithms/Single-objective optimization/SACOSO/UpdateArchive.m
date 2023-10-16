function Archive = UpdateArchive(Archive,tArchive,SwarmRBF,NDB)
% Updating database

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    % Obtain the minimum distance between swarm in SL=PSO and database
    Decs  = Archive.decs;
    D     = size(Decs,2);
    tN    = length(tArchive);
    tDecs = tArchive.decs;
    distMatrix = pdist2(Decs,SwarmRBF(:,1:D));
    distMin    = min(distMatrix,[],2);
    % Find the maximum distance in all minimum distance
    [DBPopMin,maxI]   = max(distMin);
    for k = 1 : tN
        if checkExist(Decs,tDecs(k,:))
            if length(Archive) < NDB
                Archive = [Archive,tArchive(k)];
            else
                dist = pdist2(SwarmRBF(:,1:D),tDecs(k,:));
                indPopMinK = min(dist);
                if indPopMinK < DBPopMin
                    Archive(maxI) = tArchive(k);
                    DBPopMin = indPopMinK;
                end
            end
        end
    end
end