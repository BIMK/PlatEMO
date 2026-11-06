function weightVectors = generateWeightVectors2(Nmu, epsilon)
    % Input:
    %   - Nmu: Number of weight vectors
    %   - epsilon: Small positive value (e.g., 0.01)

    % Initialize weight vectors matrix
    weightVectors = zeros(Nmu, 2);

    % Generate weight vectors
    for j = 1:Nmu
        uj1 = epsilon + (j - 1) * (1 - 2 * epsilon) / (Nmu - 1);
        uj2 = 1 - uj1;
        weightVectors(j, :) = [uj1, uj2];
    end
end
