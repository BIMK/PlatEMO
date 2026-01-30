function para = RBFCreate(ax, ay, kernel)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Jiao Liu (email: jiao.liu@ntu.edu.sg)

    [N, D] = size(ax);
    % Normlization
    xmin = min(ax, [], 1);
    xmax = max(ax, [], 1);
    ymin = min(ay, [], 1);
    ymax = max(ay, [], 1);
    ax   = 2./(repmat(xmax - xmin + eps, N, 1)) .* (ax - repmat(xmin, N, 1) + eps) - 1;
    ay   = 2./(repmat(ymax - ymin + eps, N, 1)) .* (ay - repmat(ymin, N, 1) + eps) - 1;
    r = dist(ax, ax');
    switch kernel
        case 'gaussian'
            Phi = radbas(sqrt(-log(.5))*r);
        case 'cubic'
            Phi = r.^3;
    end
    P = [ones(N, 1), ax];
    A = [Phi, P; P' zeros(D + 1, D + 1)];
    b = [ay; zeros(D + 1, size(ay, 2))];

    theta       = A \ b;
    para.alpha  = theta(1 : N, :);
    para.beta   = theta(N + 1 : end, :);
    para.xmin   = xmin;
    para.xmax   = xmax;
    para.ymin   = ymin;
    para.ymax   = ymax;
    para.nodes  = ax;   % normlization
    para.kernel = kernel;
    para.Phi    = Phi;
    para.ax     = ax;
    para.ay     = ay;
end