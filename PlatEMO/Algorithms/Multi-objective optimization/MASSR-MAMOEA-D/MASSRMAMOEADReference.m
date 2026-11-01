function reference = MASSRMAMOEADReference(Problem)
%MASSRMAMOEADReference - Construct a fixed Pareto reference set.
%
%   Request at least 10000 points and retain 10000 when more are returned.
%   Reference generation uses a fixed seed and restores the caller's
%   random-number state on exit.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    saved   = rng;
    cleanup = onCleanup(@()rng(saved));
    rng(0,'twister');
    request   = 10000;
    reference = Problem.GetOptimum(request);
    for i = 1 : 5
        if size(reference,1) >= 10000
            break;
        end
        request   = 2*request;
        reference = Problem.GetOptimum(request);
    end
    if size(reference,1) > 10000
        reference = reference(randperm(size(reference,1),10000),:);
    end
end