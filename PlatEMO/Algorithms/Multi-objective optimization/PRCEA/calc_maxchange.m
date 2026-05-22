function max_change = calc_maxchange(ideal_points,nadir_points,mean_points,gen,last_gen)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    delta_value = 1e-6 * ones(1,size(ideal_points,2));
    idx         = gen - last_gen + 1;
    rz          = abs((ideal_points(gen,:) - ideal_points(idx,:)) ./ max(abs(ideal_points(idx,:)),delta_value));
    nrz         = abs((nadir_points(gen,:) - nadir_points(idx,:)) ./ max(abs(nadir_points(idx,:)),delta_value));
    mrz         = abs((mean_points(gen,:)  - mean_points(idx,:))  ./ max(abs(mean_points(idx,:)),delta_value));
    max_change  = max([rz, nrz, mrz]);
end