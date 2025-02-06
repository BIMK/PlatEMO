function [Parent, Non0_index]= DimJud0(Parent,Problem)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Xiangyu Wang (email: xiangyu.wang@uni-bielefeld.de)

    %% feature extraction
    [~,N2] =  size(Parent);
    % density
    density = mean(Parent~=0,1);
    % calculate the statistics
    for i = 1: N2
        dim_Q75(i) = prctile(Parent(Parent(:,i)~=0,i),75);
        dim_Q50(i) = prctile(Parent(Parent(:,i)~=0,i),50);
        dim_Q25(i) = prctile(Parent(Parent(:,i)~=0,i),25);
    end
    dim_Q75(isnan(dim_Q75)) = 0;
    dim_Q50(isnan(dim_Q50)) = 0;
    dim_Q25(isnan(dim_Q25)) = 0;

    % Normalization
    ParentN = zeros(1,N2);
    % >0
    upper_index = dim_Q50>0;
    ParentN(upper_index) = dim_Q75(upper_index)./Problem.upper(upper_index);
    % <0
    lower_index = dim_Q50<0;
    ParentN(lower_index) = dim_Q25(lower_index)./Problem.lower(lower_index);
    % projection
    ParentN = sqrt(1-(ParentN-1).^2);
    % ClusteringS
    [C_indx,C_point] = ExactKmeans([ParentN',density'],4);
    [~,a] = sort(C_point(:,1));
    Is0_index1 = C_indx == a(1);
    Is0_index2 = C_indx == a(2);
    Is0_index = (Is0_index1 + Is0_index2) == 1;
    Non0_index = 1-Is0_index;
    % setting 0
    Parent(:, Is0_index) = 0;
end