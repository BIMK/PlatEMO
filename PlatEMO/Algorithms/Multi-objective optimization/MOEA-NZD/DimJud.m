function [Re, Population] = DimJud(Population, upper, lower,ReAlready)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Xiangyu Wang (email: xiangyu.wang@uni-bielefeld.de)

    [N1,N2] = size(Population);
    density = mean(Population~=0,1);
    for i = 1 : N2
        dim_Q75(i) = prctile(Population(Population(:,i)~=0,i),75);
        dim_Q50(i) = prctile(Population(Population(:,i)~=0,i),50);
        dim_Q25(i) = prctile(Population(Population(:,i)~=0,i),25);
    end
    dim_Q75(isnan(dim_Q75)) = 0;
    dim_Q50(isnan(dim_Q50)) = 0;
    dim_Q25(isnan(dim_Q25)) = 0;
    km_input = zeros(1,N2);
    upper_index = dim_Q50>0;
    km_input(upper_index) = dim_Q75(upper_index)./upper(upper_index);
    lower_index = dim_Q50<0;
    km_input(lower_index) = dim_Q25(lower_index)./lower(lower_index);
    km_input = sqrt(1-(km_input-1).^2);
    [C_indx,C_point] = ExactKmeans([km_input',density'],4);
    [~,a] = max(C_point(:,1)+C_point(:,2));
    activating_index = C_indx == a;

    %% Case 1: dim_mean > 0
    ReNow = false(1,size(ReAlready,2));
    upper_index1 = upper_index + activating_index;
    ReNow(upper_index1==2) = true;
    ReNow(ReAlready == true) = false;
    range = upper(ReNow) - dim_Q75(ReNow);
    if ~isempty(range)
        Population(:,ReNow) = rand(N1,size(range,2)).*range + repmat(dim_Q75(ReNow),N1,1);
    end
    ReAlready(ReNow ==true) = true;

    %% Case 2: dim_mean < 0
    upper_index1 =  lower_index + activating_index;
    ReNow(upper_index1==2)  = true;
    ReNow(ReAlready ==true) = false;
    range =  dim_Q25(ReNow) - lower(ReNow);
    if ~isempty(range)
        Population(:,ReNow) = -1.*rand(N1,size(range,2)).*range + repmat(dim_Q25(ReNow),N1,1);
    end
    ReAlready(ReNow ==true) = true;
    Re = ReAlready;
end