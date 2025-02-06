function weightVectors = generateWeightVectors2(Nmu, epsilon)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    % Initialize weight vectors matrix
    weightVectors = zeros(Nmu, 2);

    % Generate weight vectors
    for j = 1 : Nmu
        uj1 = epsilon + (j - 1) * (1 - 2 * epsilon) / (Nmu - 1);
        uj2 = 1 - uj1;
        weightVectors(j, :) = [uj1, uj2];
    end
end