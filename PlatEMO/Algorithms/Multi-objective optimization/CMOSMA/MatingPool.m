function [MatingPool] = MatingPool(XU,N,B)
% Construct the MatingPool

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Chao He

    MatingPool = 1 : N;
    for u = 1 :N  
        if rand < 0.9
            Q = XU(B(u,:));
        else
            Q = 1 : N;
        end
        MatingPool(u)= Q(randperm(end,1));
    end
end