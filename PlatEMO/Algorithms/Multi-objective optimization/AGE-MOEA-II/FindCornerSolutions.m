function [indexes] = FindCornerSolutions(front)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Annibale Panichella


[m,n] = size(front);

%% let's normalize the objectives
if m<=n
  indexes = 1:m;
  return
end

%% let's define the axes of the n-dimensional spaces 
W = eye(n);
[r,~]= size(W);
indexes = zeros(1,n);
for i=1:r
   [~, index] = min(Point2LineDistance(front, zeros(1,n), W(i,:)));
   indexes(i) = index;
end

end