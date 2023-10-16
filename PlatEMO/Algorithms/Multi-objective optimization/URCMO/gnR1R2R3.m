function [r1,r2,r3] = gnR1R2R3(NP1, r0)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Kangjia Qiao

    NP0 = length(r0);

    r1  = floor(rand(1, NP0) * NP1) + 1;
    for i = 1 : inf
        pos = (r1 == r0);
        if sum(pos) == 0
            break;
        else   % regenerate r1 if it is equal to r0
            r1(pos) = floor(rand(1, sum(pos)) * NP1) + 1;
        end
        if i > 1000, % this has never happened so far
            error('Can not genrate r1 in 1000 iterations');
        end
    end

    r2  = floor(rand(1, NP0) * NP1) + 1;
    for i = 1 : inf
        pos = ((r2 == r1) | (r2 == r0));
        if sum(pos) == 0
            break;
        else   % regenerate r2 if it is equal to r0 or r1
            r2(pos) = floor(rand(1, sum(pos)) * NP1) + 1;
        end
        if i > 1000, % this has never happened so far
            error('Can not genrate r2 in 1000 iterations');
        end
    end

    r3  = floor(rand(1, NP0) * NP1) + 1;
    for i = 1 : inf
        pos = ((r3 == r1) | (r3 == r0) | (r3 == r2));
        if sum(pos) == 0
            break;
        else   % regenerate r2 if it is equal to r0 or r1
            r3(pos) = floor(rand(1, sum(pos)) * NP1) + 1;
        end
        if i > 1000, % this has never happened so far
            error('Can not genrate r3 in 1000 iterations');
        end
    end
end