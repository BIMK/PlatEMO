classdef RecM2m < handle
% Crossover operator in MOEA/D-M2M

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        name    % operator name
        MaxGen  % current genetic generation
    end
    methods
        function obj = RecM2m(maxgen)
            obj.name   = 'OP2';
            obj.MaxGen = maxgen;
        end
        function OffDec = do(obj, OldChrom, r0, neighbourVector, currentGen)
            [N, D] = size(OldChrom);
            r2     = datasample(neighbourVector, 1, 'Replace', false);
            p1     = OldChrom(r0, :);
            p2     = OldChrom(r2, :);
            rc     = (2 * rand(1) - 1) * (1 - rand(1) ^ (-(1 - currentGen / (obj.MaxGen + N)) ^ 0.7));
            OffDec = p1 + rc * (p1 - p2);
        end
    end
end