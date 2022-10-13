function frontiers_picked = crowding_pick(archive, picks, mode)
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
            O = normalize([archive(index_pick).objs]);
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