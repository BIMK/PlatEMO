function Transformation = Allocation(Obj,kPoint,R)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    K     = length(R);
    [N,~] = size(Obj);
    
    % Allocation of each solution
    if K == 1
        Transformation = ones(N,1);
    else
        Transformation = zeros(N,1);
        for i = 1 : K
            T       = find(Transformation~=0);
            current = pdist2(kPoint(i,:),Obj(T,:))<=R(i);
            Transformation(T(current)) = i;
        end

        Remain = Transformation==0;
        if sum(Remain) ~= 0
            [~,transformation] = min(pdist2(Obj(Remain,:),kPoint),[],2);
            Transformation(Remain) = transformation;
        end
    end
end