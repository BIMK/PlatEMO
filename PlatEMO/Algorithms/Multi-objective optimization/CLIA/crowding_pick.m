function frontiers_picked = crowding_pick(archive, picks, mode)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

if numel(archive) < picks
    frontiers_picked = archive;
else
    O = normalize([archive.objs]);
    if strcmp(mode, 'precise')
        index_pick = true(1, numel(archive));
        while numel(find(index_pick == true)) > picks % Optimizable
            CD = NaN(1, numel(archive));
            CD(index_pick) = crowding_dist(O);
            [~, index] = min(CD);
            index_pick(index) = false;
            O = normalization([archive(index_pick).objs]);
        end
        frontiers_picked = archive(index_pick);
    elseif strcmp(mode, 'fast')
        CD = crowding_dist(O);
        [~, I] = sort(CD, 'descend');
        frontiers_picked = archive(I(1: picks));
    end
end
end

function CD = crowding_dist(O)
[N, M] = size(O);
CD = zeros(1, N);
Fmax = max(O, [], 1);
Fmin = min(O, [], 1);
for i = 1: M
    [~, Rank] = sort(O(:, i));
    CD(Rank(1)) = inf;
    CD(Rank(end)) = inf;
    for j = 2 : N - 1
        CD(Rank(j)) = CD(Rank(j)) + (O(Rank(j + 1), i) - O(Rank(j - 1), i)) / (Fmax(i) - Fmin(i));
    end
end
end