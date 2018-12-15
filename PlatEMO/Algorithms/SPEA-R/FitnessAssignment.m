function FV = FitnessAssignment(Ei,PopObj,Angle,theta)
% Calculate the local fitness value of each solution

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N = size(PopObj,1);

    %% Detect the dominance relation between each two solutions
    Dominate = false(N);
    for i = 1 : N-1
        for j = i+1 : N
            k = any(PopObj(i,:)<PopObj(j,:)) - any(PopObj(i,:)>PopObj(j,:));
            if k == 1
                Dominate(i,j) = true;
            elseif k == -1
                Dominate(j,i) = true;
            end
        end
    end
    
    %% Local strength and local raw fitness
    Sl = zeros(1,N);
    Rl = zeros(1,N);
    for i = unique(Ei)
        Local     = find(Ei==i);
        Sl(Local) = sum(Dominate(Local,Local),2);
        for j = Local
            Rl(j) = sum(Sl(Local(Dominate(Local,i))));
        end
    end
    
    %% Global strength and local raw fitness
    Sg = sum(Dominate,2);
    Rg = zeros(1,N);
    for i = 1 : N
        Rg(i) = sum(Sg(Dominate(:,i)));
    end
    
    %% Density
    D = Angle./(Angle+theta);
    
    %% Final fitness
    FV = zeros(1,N);
    for i = unique(Ei)
        Local = find(Ei==i);
        if length(Local) == 1
            FV(Local) = Rl(Local) + D(Local);
        else
            FV(Local) = Rl(Local) + D(Local) + Rg(Local);
        end
    end
end