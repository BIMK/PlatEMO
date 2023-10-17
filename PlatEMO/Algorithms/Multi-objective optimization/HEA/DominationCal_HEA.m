function [non_dominated, hd] = DominationCal_HEA(Objs, zmin, zmax, T)
% The calculation of hyper-dominance degree

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhe Liu

    %% Normalization the solutions
    [Objs_num, M] = size(Objs);
    objs_dominated = zeros(1, Objs_num);
    hd = zeros(1, Objs_num);
    NormalizedObj = Objs ./ (zmax - zmin);
    
    %% Calculate hyper-dominance degree
    for i = 1: Objs_num
        err = NormalizedObj(i, :) - NormalizedObj;
        eq = zeros(Objs_num, 1);
        max_err = max(err, [],  2);
        min_err = min(err, [],  2);
        H = -min_err ./ max_err;
        for j = i + 1: Objs_num
            for k = 1: M
                if err(j, k) ~= 0
                    break
                end
                if k == M
                    eq(j) = 1;
                end
            end
        end
        objs_dominated(eq == 1) = 1;
        H(eq == 1) = -inf;  
        H(max_err <= 0) = inf;
        hd(i) = min(H);
    end
    
    %% Identify dominated solutions
    objs_dominated(hd < T) = 1;
    non_dominated = 1 - objs_dominated; 
end