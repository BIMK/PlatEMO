function [X,Bool] = Compare(Pop,priority,current_cons,constraint_handing)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    if current_cons == 0
        CV = zeros(2,1);
    else
        CV = Pop.cons;
        CV = CV(:,priority(1:current_cons));
        CV = sum(max(0,CV),2);
    end
    if constraint_handing == 0
        if CV(1) < CV(2)
            Bool = true;
            X    = Pop(1);
        elseif CV(1) > CV(2)
            Bool = false;
            X    = Pop(2);
        else
            k = any(Pop(1).objs<Pop(2).objs) - any(Pop(1).objs>Pop(2).objs);
            if k == 1 
                Bool = true;
                X = Pop(1);
            elseif k == -1 || k == 0
                Bool = false;
                X    = Pop(2);
            end
        end
    end
end