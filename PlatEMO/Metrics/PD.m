function score = PD(Population,~)
% <max> <multi/many> <real/integer/label/binary/permutation> <large/none> <constrained/none> <expensive/none> <multimodal/none> <sparse/none> <dynamic/none>
% Pure diversity

%------------------------------- Reference --------------------------------
% H. Wang, Y. Jin, and X. Yao, Diversity assessment in many-objective
% optimization, IEEE Transactions on Cybernetics, 2017, 47(6): 1510-1522.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    PopObj = Population.best.objs;
    if isempty(PopObj)
        score = nan;
    else
        C = false(length(Population));
        C(logical(eye(size(C)))) = true;
        D = pdist2(Population.objs,Population.objs,'minkowski',0.1);
        D(logical(eye(size(D)))) = inf;
        score = 0;
        for k = 1 : length(Population)-1
            while true
                [d,J] = min(D,[],2);
                [~,i] = max(d);
                if D(J(i),i) ~= -inf
                    D(J(i),i) = inf;
                end
                if D(i,J(i)) ~= -inf
                    D(i,J(i)) = inf;
                end
                P = any(C(i,:),1);
                while ~P(J(i))
                    newP = any(C(P,:),1);
                    if P == newP
                        break;
                    else
                        P = newP;
                    end
                end
                if ~P(J(i))
                    break;
                end
            end
            C(i,J(i)) = true;
            C(J(i),i) = true;
            D(i,:)    = -inf;
            score     = score + d(i);
        end
    end
end