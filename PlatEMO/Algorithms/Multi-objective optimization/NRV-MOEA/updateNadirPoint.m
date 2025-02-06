function [Znadir,extremPoint] = updateNadirPoint(PopObj,Zideal,Znadir,extremPoint)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    for i = 1 : size(PopObj,1)
        for j = 1 : size(PopObj,2)
            asf1 = asfFunction(PopObj(i,:),j,Zideal,Znadir);
            asf2 = asfFunction(extremPoint(j,:),j,Zideal,Znadir);
            if asf1 < asf2
                extremPoint(j,:) = PopObj(i,:);
            end
        end
    end

    % update the nadir point
    M    = size(PopObj,2);
    temp = extremPoint - repmat(Zideal,M,1);

    if rank(temp) == size(temp,1)
        u  = ones(M,M);
        al = inv(temp)*u;
        for j = 1 : M
            aj = 1./al(j,1) + Zideal(j);
            if (aj > Zideal(j)) && ~isinf(aj) && ~isnan(aj)
                Znadir(j) = aj;
            else
                break;
            end
        end
    else
        zmax = max(PopObj,[],1);
        Znadir = zmax;
    end
end

function maxValue = asfFunction(sol,index,Zideal,Znadir)
    epsilon  = 1.0E-6;
    maxValue = 0;
    for i = 1 : size(sol,2)
        val = abs((sol-Zideal)./(Znadir - Zideal));
        if index ~= i
            val = val/epsilon;
        end
        if val > maxValue
            maxValue = val;
        end
    end
end