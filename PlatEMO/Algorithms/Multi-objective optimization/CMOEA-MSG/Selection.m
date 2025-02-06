function [X,oc] = Selection(Population,current_cons,Fitness,priority)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    cons = Population.cons;
    cons = cons(:,priority);
    if current_cons == 0
        oc = 0;
        [~,temp] = sort(Fitness);
        X = [Population(temp(floor(end/10))),Population(randi(100,1,2))];
    elseif current_cons == 1
        temp = find(cons(:,current_cons)>0);
        if isempty(temp)
            oc = 0;
            [~,temp] = sort(Fitness);
            X = Population(temp(1));
        else
            oc = 1;
            [~,I] = sort(Fitness(temp));
            X = Population(temp(I(end)));
        end
    else
        temp  = find(max(cons(:,1:current_cons-1),[],2)<0);
        index = find(cons(temp,current_cons)>0);
        if isempty(index)
            oc = 0;
            [~,temp] = sort(Fitness);
            X = Population(temp(1));
        else
            oc    = 1;
            temp  = temp(index);
            [~,I] = sort(Fitness(temp));
            X     = Population(temp(I(end)));
        end
    end
end