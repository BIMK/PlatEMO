function [MergeDec] = MergeSubpop(eachDec,Index,cycle,NumEsp,mu,Problem)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yingwei Li

    MergeDec = zeros(mu,Problem.D);
    for j = 1 : cycle
        for i = 1 : NumEsp  
            MergeDec(:,Index{j}==i) = eachDec{i}; 
        end
    end
end