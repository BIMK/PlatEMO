function Level = LevelSort(Population,n)
% Ssort the level of Population

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming

    Zmax     = max(Population.objs,[],1);
    Zmin     = min(Population.objs,[],1);
    interval = (Zmax-Zmin)./n;
    Level    = zeros(length(Population),1);
    objs     = Population.objs;
    for i = 1 : length(Population)
        t = 0;
        leveled = 0;
        obj = objs(i,:);
        while leveled == 0
            t = t+1;
            leveled = 1;
            for j = 1 : size(objs,2)
                if obj(1,j)>Zmin(j)+(t+1)*interval(1,j)
                    leveled = 0;
                    break;
                end
            end
        end
        Level(i) = t;
    end 
end