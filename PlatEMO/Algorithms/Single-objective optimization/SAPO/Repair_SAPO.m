function x = Repair_SAPO(x, xL, xU)
%% Repair the solution that violates the boundary constraints (reflection)

    for i = 1 : size(x, 1)
        indexLower = find(x(i, :) < xL);
        indexUper = find(x(i, :) > xU);
        x(i, indexLower) = min(xU(indexLower), 2 * xL(indexLower) - x(i, indexLower));
        x(i, indexUper) = max(xL(indexUper), 2 * xU(indexUper) - x(i, indexUper));
    end
end