function Archive = UpdateArchive(Archive,N)
% Update the offline archive in PICEA-g

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Archive = Archive(NDSort(Archive.objs,1)==1);
    Archive = Truncation(Archive,N);
end

function Population = Truncation(Population,N)
% Select part of the solutions by truncation

    Distance = pdist2(Population.objs,Population.objs);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,length(Population));
    while sum(Del) < length(Population)-N
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
    Population(Del) = [];
end