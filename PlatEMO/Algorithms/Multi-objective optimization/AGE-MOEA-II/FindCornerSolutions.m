function [indexes] = FindCornerSolutions(front)
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