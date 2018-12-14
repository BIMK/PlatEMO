function delta = getDelataMinfor(Population,ObjectiveSet1,ObjectiveSet2)
% Calculate the objective subset with KMOSS greedy method

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Liangli Zhen

   delta   = 0;
   PopSize = size(Population.decs,1);
   PopObj  = Population.objs;
   
   %% Calculate the weakly dominated solution pairs
   DominateRelation = zeros(PopSize,PopSize);
   for i = 1 : PopSize
       DominateRelation(:, i) = ~any(repmat(PopObj(i,ObjectiveSet1),PopSize,1)>PopObj(:,ObjectiveSet1),2);
   end
   [row, col] = find(DominateRelation==1);
   
   %% Compute the minimum value that makes the dominance relation unchanged
   delta = max(delta,max(max(PopObj(col,ObjectiveSet2)-PopObj(row,ObjectiveSet2))));
end