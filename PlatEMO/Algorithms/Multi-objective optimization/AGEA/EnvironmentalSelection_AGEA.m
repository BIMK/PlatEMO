function [Grid, GridIndex] = EnvironmentalSelection_AGEA(Grid, zmin, gmax, div)
PopObj = Grid.objs;
[N, M] = size(PopObj);
GridIndex = zeros(N, M);
GridCorner = zeros(N, M);
delta = ones(N, M);    % The blow-up factor
side_length    = (gmax - zmin) / (div - 1);
lower_boundary    = zmin - 0.5 * side_length;
for i = 1: N
    for j = 1: M
        GridIndex(i, j) = floor((PopObj(i, j) - lower_boundary(j)) / side_length(j));
        GridCorner(i, j) = lower_boundary(j) + side_length(j) *  GridIndex(i, j);
        if GridIndex(i, j) == 0
            delta(i, j) = 1e+6;
        end
    end
end    
normalized_obj = (PopObj - zmin) ./ (gmax - zmin);
normalized_GridCorner = (GridCorner - zmin) ./ (gmax - zmin);
Loser = false(N, 1);
fitness = sum(((normalized_obj - normalized_GridCorner) .* delta) .^ 2, 2);
for i = 1: N - 1
    for j = i + 1: N
        if GridIndex(i,:) == GridIndex(j,:)
            if fitness(i) <= fitness(j)
                Loser(j) = true;
            else
                Loser(i) = true;
            end
        end
    end
end
Grid = Grid(~Loser);
GridIndex = GridIndex(~Loser, :);
end