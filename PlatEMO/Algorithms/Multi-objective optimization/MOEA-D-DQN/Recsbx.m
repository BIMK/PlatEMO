classdef Recsbx
% Simulated binary crossover

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        name        % operator name
        XOVR        % cross probability
        Half        % half uniform cross
        n           % population distribution
    end
    
    methods
        function obj = Recsbx()
            obj.name = 'OP1';
            obj.XOVR = 0.7;
            obj.Half = true;
            obj.n    = 20;
        end
        function off = do(obj, OldChrom, r0, neighbourVector)
            r1   = datasample(neighbourVector, 1);
            p0   = OldChrom(r0, :);
            p1   = OldChrom(r1, :);
            D    = size(p1, 2);
            mu   = rand(1, D);
            beta = zeros(1, D);
            idx  = mu <= 0.5;
            beta(idx)  = (2 * mu(idx)).^(1 / (obj.n + 1));
            beta(~idx) = (2 - 2 * mu(~idx)).^(-1 / (obj.n + 1));
            beta = beta .* randsample([-1, 1], 1, true, [0.5, 0.5]);
            idx  = rand(1, D) < 0.5;
            beta(idx) = 1;
            if rand < 0.5
                off = (p0 + p1) / 2 + beta .* (p0 - p1) / 2;
            else
                off = (p0 + p1) / 2 - beta .* (p0 - p1) / 2;
            end
        end
    end
end