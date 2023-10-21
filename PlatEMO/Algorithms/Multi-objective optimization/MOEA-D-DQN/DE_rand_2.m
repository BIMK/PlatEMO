classdef DE_rand_2
% Differential evolution type 2

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
        F       % scaling factor
    end
    methods
        function obj = DE_rand_2()
            obj.name = 'OP4';
            obj.F    = 0.5;
        end
        function v = do(obj, OldChrom, r0, neighbourVector)
            x  = OldChrom(r0, :);
            r1 = neighbourVector(1);
            r2 = neighbourVector(2);
            r3 = neighbourVector(3);
            r4 = neighbourVector(4);
            v  = x + obj.F * (OldChrom(r1, :) - OldChrom(r2, :) + OldChrom(r3, :) - OldChrom(r4, :));
        end
    end
end