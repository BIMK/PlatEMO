function [x_pca, m] = PCAMyself(x)
% Feature extraction by PCA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [COEFF, SCORE, latent]=pca(x);
    answer=0;s1=0;s2=sum(latent);count=1;
    while answer < 0.95
        s1=s1+latent(count);
        answer=s1/s2;
        count=count+1;
    end
    m=COEFF(:,1:count-1);
    x_pca=x*m;
end