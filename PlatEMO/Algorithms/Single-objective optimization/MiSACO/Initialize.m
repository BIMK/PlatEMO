function [v_r, v_o, v_c, f, rank_out] = Initialize(Problem)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Jiao Liu (email: jiao.liu@ntu.edu.sg)

    up_r  = Par.up_r;
    dn_r  = Par.dn_r;
    up_o  = Par.up_o;
    dn_o  = Par.dn_o;
    l     = Par.l;
    k     = Par.k;
    len_r = Par.len_r;
    len_o = Par.len_o;
    len_c = Par.len_c;
    func_num = Par.func_num;
    shift    = Par.shift;
    
    v_r = repmat((up_r - dn_r), k, 1) .* lhsdesign(k, len_r) + repmat(dn_r, k, 1);
    v_o = repmat((up_o - dn_o), k, 1) .* lhsdesign(k, len_o) + repmat(dn_o, k, 1);
    v_c = ceil(repmat(l, k, 1) .* lhsdesign(k, len_c));
    f   = Fitness_func(v_r, v_o, v_c,func_num,shift);
    
    [~,rank_v]   = sort(f);
    [~,rank_out] = sort(rank_v);
end