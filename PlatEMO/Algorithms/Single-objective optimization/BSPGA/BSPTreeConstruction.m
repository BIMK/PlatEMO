function [Population,T,best] = BSPTreeConstruction(Problem,P,T,best)
% Construction of the binary space partition tree

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Population = [];
    for i = 1 : size(P,1)
        if isempty(T)
            Population = [Population,Problem.Evaluation(P(i,:))];
            T    = NODE(P(i,:),1);
            best = struct('dec',P(i,:),'fitness',FitnessSingle(Population(end)),'level',1);
        else
            SubT = T;
            while ~isempty(SubT.left)
                if P(i,SubT.level) == 0
                    SubT = SubT.left;
                else
                    SubT = SubT.right;
                end
            end
            if SubT.level <= size(P,2)
                if P(i,SubT.level) == SubT.node(SubT.level)
                    if sum(P(i,:)~=best.dec) < sum(SubT.node~=best.dec)
                        P(i,SubT.level) = ~P(i,SubT.level);
                    else
                        SubT.node(SubT.level) = ~SubT.node(SubT.level);
                        Population = [Population,Problem.Evaluation(SubT.node)];
                        if FitnessSingle(Population(end)) < best.fitness
                            best = struct('dec',Population(end).dec,'fitness',FitnessSingle(Population(end)),'level',SubT.level+1);
                        end
                    end
                end
                if SubT.node(SubT.level) == 0
                    SubT.left  = NODE(SubT.node,SubT.level+1);
                    SubT.right = NODE(P(i,:),SubT.level+1);
                else
                    SubT.left  = NODE(P(i,:),SubT.level+1);
                    SubT.right = NODE(SubT.node,SubT.level+1);
                end
                Population = [Population,Problem.Evaluation(P(i,:))];
                if FitnessSingle(Population(end)) < best.fitness
                    best = struct('dec',Population(end).dec,'fitness',FitnessSingle(Population(end)),'level',SubT.level+1);
                end
            end
        end
    end
end