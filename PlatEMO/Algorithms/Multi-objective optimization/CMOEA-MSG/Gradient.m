function Gradient = Gradient(Problem,Pop,current_cons,priority,oc)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    if oc == 1
        % Gradient of a constraint
        [~,ConGrad] = Problem.CalGrad(Pop.dec);
        Gradient    = ConGrad(priority(current_cons),:);
    else
        % Gradient of all objectives
        ObjGrad  = Problem.CalGrad(Pop.dec);
        Gradient = sum(ObjGrad,1);
        conflict = any(ObjGrad<0,1) & any(ObjGrad>0,1);
        Gradient(conflict) = 0;
    end
end