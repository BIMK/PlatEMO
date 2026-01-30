function y = RBFInterp(x, para)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Jiao Liu (email: jiao.liu@ntu.edu.sg)

    ax = para.nodes;
    nx = size(x, 1);
    
    xmin = para.xmin;
    xmax = para.xmax;
    ymin = para.ymin;
    ymax = para.ymax;
    % normalization
    x = 2./(repmat(xmax - xmin + eps, nx, 1)) .* (x - repmat(xmin, nx, 1) + eps) - 1;
    r = dist(x, ax');
    switch para.kernel
        case 'gaussian'
            Phi = radbas(sqrt(-log(.5))*r);
        case 'cubic'
            Phi = r.^3;
    end
    y = Phi * para.alpha + [ones(nx, 1), x] * para.beta;
    % renormalization
    y = repmat(ymax - ymin + eps, nx, 1)./2 .* (y + 1) + repmat(ymin - eps, nx, 1);
end