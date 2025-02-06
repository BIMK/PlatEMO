function [priority,Feasible_rate] = Constraint_priority(Population,priority,current_cons)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
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
    if isempty(priority)
        [~,temp] = sort(Feasible_rate);
        priority = temp;
    elseif current_cons+1 < size(Population.cons,2)
        Feasible_rate = Feasible_rate(:,priority);
        [~,temp] = sort(Feasible_rate(current_cons+1:end));
        priority(current_cons+1:end) = priority(temp +current_cons)   ;
    end
end