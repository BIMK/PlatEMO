function TransferredDec = domainAdaption(TDecs, ODecs, x)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    dim     = size(ODecs, 2);     
    coeff_t = pca(TDecs, 'NumComponents', floor(dim * 0.5));
    coeff_o = pca(ODecs, 'NumComponents', floor(dim * 0.5));

    orth_coeff_t = orth(coeff_t);
    orth_coeff_o = orth(coeff_o);

    Xb = orth_coeff_t * orth_coeff_o' * orth_coeff_t;

    o_pca_sa    = ODecs * Xb;
    o_pca_sa_re = o_pca_sa * coeff_t';
    max_v       = max(max(o_pca_sa_re));
    min_v       = min(min(o_pca_sa_re));

    TransferredDec = (o_pca_sa_re(x, :) - min_v) ./ (max_v - min_v);
end