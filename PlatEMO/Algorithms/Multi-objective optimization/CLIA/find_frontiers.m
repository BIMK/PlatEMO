function [F_index, NF_index] = find_frontiers(PopObj, PopCon)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

	perm_index = randperm(size(PopObj, 2));
	PopObj = PopObj(:, perm_index);
	Infeasible = any(PopCon>0, 2);
	PopObj(Infeasible,:) = repmat(max(PopObj,[],1),sum(Infeasible),1) + repmat(sum(max(0,PopCon(Infeasible,:)),2),1,size(PopObj,2));
 [PopObj, ~, Loc] = unique(PopObj,'rows');
 [N,M] = size(PopObj);
 [PopObj,rank] = sortrows(PopObj);
 FrontNo = inf(1,N);
 MaxFNo = 1;
 for i = 1 : N
 if FrontNo(i) == inf
 Dominated = false;
 for j = i-1 : -1 : 1
 if FrontNo(j) == MaxFNo
 m = 2;
 while m <= M && PopObj(i, m) >= PopObj(j, m)
 m = m + 1;
 end
 Dominated = m > M;
 if Dominated || M == 2
 break;
 end
 end
 end
 if ~Dominated
 FrontNo(i) = MaxFNo;
 end
 end
 end
 FrontNo(rank) = FrontNo;
 FrontNo = FrontNo(Loc);
 F_index = find(FrontNo == 1);
 NF_index = setdiff(1: size(PopObj, 1), F_index);
end