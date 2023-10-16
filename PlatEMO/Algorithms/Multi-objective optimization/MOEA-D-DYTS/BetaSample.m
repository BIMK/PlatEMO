function op = BetaSample(pds)
% Thompson Sampling-based operator selection

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    n = size(pds,2);
    ops = ones(1,n);
    for i = 1:size(pds,2)
        ops(i) = random(pds(i));
    end
    [~,op] = max(ops);
end