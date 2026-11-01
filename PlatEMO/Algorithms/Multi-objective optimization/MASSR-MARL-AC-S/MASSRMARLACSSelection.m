function [Population,Dec,Mask,FitnessSpea2] = MASSRMARLACSSelection(Population,Dec,Mask,N)
%MASSRMARLACSSelection - Select sparse solutions with SPEA2 fitness.
%
%   Remove duplicate objective vectors, then apply strength-and-density
%   selection and distance truncation. Population, latent Dec, Mask, and
%   returned fitness remain aligned.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    % Remove duplicate objectives while retaining aligned decisions
    [~,uni]    = unique(Population.objs,'rows');
    Population = Population(uni);
    Dec        = Dec(uni,:);
    Mask       = Mask(uni,:);
    N          = min(N,length(Population));

    Fitness = MASSRMARLACSFitness(Population.objs);

    % Select by SPEA2 fitness and truncate crowded survivors
    Next = Fitness < 1;
    if sum(Next) < N
        [~,Rank]        = sort(Fitness);
        Next(Rank(1:N)) = true;
    elseif sum(Next) > N
        Del             = Truncation(Population(Next).objs,sum(Next)-N);
        Temp            = find(Next);
        Next(Temp(Del)) = false;
    end
    Population   = Population(Next);
    Fitness      = Fitness(Next);
    Dec          = Dec(Next,:);
    Mask         = Mask(Next,:);
    FitnessSpea2 = Fitness;
end

function Del = Truncation(PopObj,K)
%Truncation - Identify crowded solutions for removal.
%
%   Repeatedly delete the solution with the lexicographically smallest
%   sorted distance vector until K deletions have been selected.

    Distance                                 = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del                                      = false(1,size(PopObj,1));
    while sum(Del) < K
        Remain               = find(~Del);
        Temp                 = sort(Distance(Remain,Remain),2);
        [~,Rank]             = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end
