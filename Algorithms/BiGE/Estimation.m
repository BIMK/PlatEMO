function BiObj = Estimation(PopObj,r)
% Estimate the proximity and crowding degree of each solution

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N = size(PopObj,1);
    
    %% Proximity estimation
    fmax   = repmat(max(PopObj,[],1),N,1);
    fmin   = repmat(min(PopObj,[],1),N,1);
    PopObj = (PopObj-fmin)./(fmax-fmin);
    fpr    = sum(PopObj,2);
    
    %% Crowding degree estimation
    d     = pdist2(PopObj,PopObj);
    d(logical(eye(length(d)))) = inf;
    fprm  = repmat(fpr,1,N);
    case1 = d<r & fprm<=fprm';
    case2 = d<r & fprm>fprm';
    sh        = zeros(N);
    sh(case1) = (0.5*(1-d(case1)/r)).^2;
    sh(case2) = (1.5*(1-d(case2)/r)).^2;
    fcd   = sqrt(sum(sh,2));
    BiObj = [fpr,fcd];
end