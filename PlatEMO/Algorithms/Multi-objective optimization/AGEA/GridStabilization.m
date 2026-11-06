function gmax = GridStabilization(NDPop, zmin, gmax, div)
    [~, M] = size(NDPop.objs);
    zmax = max(max(NDPop.objs, [], 1), zmin + 1e-10);
    side_length = (gmax - zmin) / (div - 1);
    for i = 1: M
        if abs(zmax(i) - gmax(i)) > 0.5*side_length(i)
            gmax(i) = zmax(i);
        end
    end       
end