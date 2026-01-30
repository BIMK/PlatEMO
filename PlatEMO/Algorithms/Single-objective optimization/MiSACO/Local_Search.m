function out_best = Local_Search(best_Arc_r,best_Arc_f,local_up,local_dn,best_r)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Jiao Liu (email: jiao.liu@ntu.edu.sg)

    global rbfnet_best 
    
    if any(isnan(best_r))
        best_r = local_dn + (local_up - local_dn) * rand(1, length(local_up));
    end
    best_r = max([best_r;local_dn]);
    best_r = min([best_r;local_up]);
    
    [~,distinct1] = unique(round(best_Arc_r*1e10)/1e10,'rows');
    [~,distinct2] = unique(round(best_Arc_r(:,1)*1e10)/1e10,'rows');
    distinct      = intersect(distinct1,distinct2);
    best_Arc_r    = best_Arc_r(distinct,:);
    best_Arc_f    = best_Arc_f(distinct,:);
    
    rbfnet_best = RBFCreate(best_Arc_r, best_Arc_f(:,1),'cubic');
    
    options = optimset('Algorithm', 'interior-point', 'MaxFunEvals', 300, 'Display', 'off');
    obj(best_r');
    out_best = (fmincon(@obj, (best_r)', [], [], [], [], local_dn, local_up, [], options))';
end

function y = obj(x)
    global rbfnet_best
    y = RBFInterp(x',rbfnet_best);
end