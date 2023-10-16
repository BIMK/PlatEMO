function [priority,Feasible_rate] = Constraint_priority(Population)
% Determine constraint-handling priority

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Feasible_rate = zeros(1,size(Population.cons,2));
    for j = 1 : size(Population.cons,2)
        CV = Population.cons;
        CV = CV(:,j);
        Feasible_rate(1,j) = length(find(CV<=0))/size(Population,2);
    end
    [~,temp] = sort(Feasible_rate);
    priority = temp;          
end