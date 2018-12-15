function Output = GetOutput(PopObj,RefPoint)
% Character with two types of solutions, 0 or 1

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He

    N = size(PopObj,1);
    Output = true(N,1);
    for i = 1 : size(RefPoint,1)
        Output = Output & any(PopObj<=repmat(RefPoint(i,:),N,1),2);
    end
end