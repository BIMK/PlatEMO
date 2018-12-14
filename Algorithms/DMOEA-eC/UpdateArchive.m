function [Archive,znad] = UpdateArchive(Archive,N)
% Update the archive in DMOEA-¦ÅC

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    %% Detect the non-dominated solutions
    Archive = Archive(NDSort(Archive.objs,1)==1);
    AObj    = Archive.objs;
    
    %% Select the extreme solutions
    Choose           = false(1,length(Archive)); 
    [~,extreme1]     = max(AObj,[],1);
    [~,extreme2]     = min(AObj,[],1);
    Choose(extreme1) = true;
    Choose(extreme2) = true;
    
    %% Select other solutions by truncation
    if sum(Choose) > N
        selected = find(Choose);
        Choose   = selected(randperm(length(selected),N));
    else
        Distance = pdist2(AObj,AObj);
        Distance(logical(eye(length(Distance)))) = inf;
        while sum(Choose) < N && ~all(Choose)
            unSelected = find(~Choose);
            [~,x]      = max(min(Distance(~Choose,Choose),[],2));
            Choose(unSelected(x)) = true;
        end
    end
    Archive = Archive(Choose);
    
    %% Update the nadir point
    znad = max(Archive.objs,[],1);
end