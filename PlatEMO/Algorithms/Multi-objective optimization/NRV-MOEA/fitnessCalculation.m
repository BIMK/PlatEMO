function [fitness] = fitnessCalculation(PopObjn,Global,Hyperplane_bp)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [mapPop,Hyperplane] = vertmap(PopObjn,PopObjn,Hyperplane_bp,Global);
    T = clusterdata(mapPop,'maxclust',size(PopObjn,1),'distance','euclidean','linkage','ward');
    for c = 1 : size(PopObjn,1)
        current = find(T == c);
        pn  = length(current);
        Ref = sum(mapPop(current,:),1)/pn;
        d12 = zeros(pn,1);
        for pc = 1 : pn
            d1 = norm(Ref-mapPop(current(pc),:));
            d2 = -(PopObjn(current(pc),:)*Hyperplane-1)./sqrt(sum(Hyperplane.^2));
            d12(pc) = d1-d2;
        end
        [~,ct] = min(d12);
        choose = current(ct);
        fitness(choose) = d12;
        fitness(choose) = sum(PopObjn(choose,:).^2);
    end
end