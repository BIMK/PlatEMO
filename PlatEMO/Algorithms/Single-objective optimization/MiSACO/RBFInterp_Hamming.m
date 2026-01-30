function [y, sigma2] = RBFInterp_Hamming(x,len_c, up, dn, para)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Jiao Liu (email: jiao.liu@ntu.edu.sg)

    [~,D] = size(x);
    ax    = para.nodes;
    nx    = size(x, 1);  
    
    ymin  = para.ymin;
    ymax  = para.ymax;
    ax_nc = ax(:,1:D-len_c);
    ax_c  = ax(:,D-len_c+1:D);
    x_nc  = x(:,1:D-len_c);
    x_nc  = (x_nc - dn)./(up - dn) - 0.5; 
    x_c   = x(:,D-len_c+1:D);
    
    dis_euc   = sum(abs(x_nc - ax_nc),2)';
    dis_gower = sum( ax_c~= x_c ,2)';
    
    r = dis_gower + dis_euc;
    switch para.kernel
        case 'gaussian'
            Phi = radbas(sqrt(-log(.5))*r);
        case 'cubic'
            Phi = r.^3;
    end
    
    y = Phi * para.alpha;
    % renormalization
    y = repmat(ymax - ymin, nx, 1)./2 .* (y + 1) + repmat(ymin, nx, 1);
    
    switch para.kernel
        case 'gaussian'
            sigma2 = 1 / sqrt(2 * pi) - Phi*(para.Phi\(Phi'));
        case 'cubic'
            sigma2 = - Phi*(para.Phi\(Phi'));
    end
    sigma2 = abs(diag(sigma2));
end