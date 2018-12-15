function objective_subset = kEMOSS(Population,k)
% Calculate the objective subset with greedy k-EMOSS method

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Liangli Zhen

   objective_set   = 1:size(Population.objs,2);
   selected_subset = [];
   while length(selected_subset) < k
       unselected_subset = setdiff(objective_set,selected_subset);
       errors = zeros(length(unselected_subset),1);
       for i = 1:length(unselected_subset)
           errors(i,1) = getDelataMinfor(Population, [selected_subset, unselected_subset(i)], objective_set);
       end
       [~, v] = min(errors);
       selected_subset = [selected_subset unselected_subset(v(1))];
   end
   objective_subset = selected_subset;
end