function [Population, Population_hd] = EnvironmentalSelection_HEA(Solutions, hd, zmin, zmax, N, W)
PopObj = Solutions.objs;
[P, ~] = size(PopObj);
zmin       = min(zmin, min(PopObj,[],1));
nf         = zmax - zmin;  % normalization factor
nf(nf < 1e-6)   = 1e-6;                     % Prevent dividing by zero
PopObj       = (PopObj - zmin) ./ nf;
Fitness = pdist2(PopObj,W,'cosine');
[~, index] = min(Fitness, [], 1);
non_dominated = zeros(P, 1);
for i = 1: N
    non_dominated(index(i)) = 1;
end
Population = Solutions(non_dominated == 1);
Population_hd = hd(non_dominated == 1);
K = N - sum(non_dominated);
for i = 1: K
    [~, Maxhd_index] = max(hd);
    Population = [Population, Solutions(Maxhd_index)];
    Population_hd = [Population_hd, hd(Maxhd_index)];
    hd(Maxhd_index) = -inf;
end
end