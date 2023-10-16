function [non_dominated, hd] = DominationCal_HEA(Objs, zmin, zmax, T)
[Objs_num, M] = size(Objs);
objs_dominated = zeros(1, Objs_num);
hd = zeros(1, Objs_num);
NormalizedObj = Objs ./ (zmax - zmin);
for i = 1: Objs_num
    err = NormalizedObj(i, :) - NormalizedObj;
    eq = zeros(Objs_num, 1);
    max_err = max(err, [],  2);
    min_err = min(err, [],  2);
    H = -min_err ./ max_err;
    for j = i + 1: Objs_num
        for k = 1: M
            if err(j, k) ~= 0
                break
            end
            if k == M
                eq(j) = 1;
            end
        end
    end
    objs_dominated(eq == 1) = 1;
    H(eq == 1) = -inf;  
    H(max_err <= 0) = inf;
    hd(i) = min(H);
end
objs_dominated(hd < T) = 1;
non_dominated = 1 - objs_dominated; 
end