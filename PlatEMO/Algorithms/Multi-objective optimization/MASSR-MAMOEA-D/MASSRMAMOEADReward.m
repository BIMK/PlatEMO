function [reward,tracker] = MASSRMAMOEADReward(value,tracker)
%MASSRMAMOEADReward - Compute the best-so-far IGD triangle reward.
%
%   An empty tracker initializes the reference IGD and progress. A new best
%   IGD receives half the increase in squared normalized progress; all other
%   steps receive zero.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    reward = 0;
    if isempty(tracker)
        tracker = struct('initial',max(value,1e-12),'best',value,'progress',0);
    elseif value < tracker.best
        progress         = (tracker.initial-value)/tracker.initial;
        reward           = 0.5*(progress^2-tracker.progress^2);
        tracker.best     = value;
        tracker.progress = progress;
    end
end