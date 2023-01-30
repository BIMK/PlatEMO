function gBest = gBest_get(AA)
% The gBest (global best) updating strategy of S-ECSO

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [M,~]     = size(AA.decs);
    r         = rand(1,size(AA.objs,2));
    r_matr    = repmat(r,M,1);
    f_gBest   = sum(r_matr.*AA.objs,2)/sum(r);
    [~,index] = min(f_gBest);
    B         = AA.decs;
    gBest     = B(index,:);
end