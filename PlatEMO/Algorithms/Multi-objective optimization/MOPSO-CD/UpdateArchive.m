function Archive = UpdateArchive(Archive,N)
% Update the archive

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Find the non-dominated solutions
    Archive = Archive(NDSort(Archive.objs,1)==1);
    
    %% Truncate the archive according to the crowding distances
    while length(Archive) > N
        [~,rank] = sort(CrowdingDistance(Archive.objs));
        Archive(rank(randi(ceil(length(rank)*0.1)))) = [];
    end
    [~,rank] = sort(CrowdingDistance(Archive.objs),'descend');
    Archive  = Archive(rank);
end