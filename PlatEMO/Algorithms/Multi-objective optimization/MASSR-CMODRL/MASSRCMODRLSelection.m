function [Population,Fitness,Next] = MASSRCMODRLSelection(Population,N,mode,epsilon)
%MASSRCMODRLSelection - Perform SPEA2 selection with constraint modes.
%
%   mode 0 ignores constraints, mode 1 uses constraint dominance, and mode 2
%   treats violations up to epsilon as feasible. Next identifies survivors
%   in the input candidate set.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Obj   = Population.objs;
    count = size(Obj,1);

    % Construct the violation values used by the selected mode
    CV = zeros(count,1);
    if mode > 0
        CV = sum(max(0,Population.cons),2);
    end
    if mode == 2
        CV(CV<=epsilon) = 0;
    end

    % Strength and raw fitness under the selected dominance rule
    dominates = false(count);
    for i = 1 : count-1
        for j = i+1 : count
            if CV(i) < CV(j)
                dominates(i,j) = true;
            elseif CV(i) > CV(j)
                dominates(j,i) = true;
            elseif all(Obj(i,:)<=Obj(j,:)) && any(Obj(i,:)<Obj(j,:))
                dominates(i,j) = true;
            elseif all(Obj(j,:)<=Obj(i,:)) && any(Obj(j,:)<Obj(i,:))
                dominates(j,i) = true;
            end
        end
    end
    strength = sum(dominates,2);
    raw      = dominates'*strength;

    % Objective-space density and survivor truncation
    distance                = MASSRFeatures.distance(Obj,Obj);
    distance(1:count+1:end) = inf;
    ordered                 = sort(distance,2);
    Fitness                 = (raw+1./(ordered(:,max(1,floor(sqrt(count))))+2))';
    N                       = min(N,count);
    Next                    = Fitness < 1;
    if sum(Next) < N
        [~,order]        = sort(Fitness);
        Next(order(1:N)) = true;
    elseif sum(Next) > N
        while sum(Next) > N
            keep                 = find(Next);
            local                = sort(distance(keep,keep),2);
            [~,order]            = sortrows(local);
            Next(keep(order(1))) = false;
        end
    end
    Population = Population(Next);
    Fitness    = Fitness(Next);
end