function [XU] = Associate(Population,W,N)
% Associate each solution with a neuron

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Chao He

     A  = 1 : N;
     U  = 1 : N;
     XU = zeros(1,N);
     for i = 1 : N
         x = randi(length(A));
         [~,u] = min(pdist2(Population(A(x)).dec,W(U,:)));
         XU(U(u)) = A(x);
         A(x)     = [];
         U(u)     = [];
     end
end