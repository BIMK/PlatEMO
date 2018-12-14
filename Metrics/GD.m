function Score = GD(PopObj,PF)
% <metric> <min>
% Generational distance

%------------------------------- Reference --------------------------------
% D. A. Van Veldhuizen, Multiobjective evolutionary algorithms:
% Classifications, analyses, and new innovations, Ph.D. thesis, Department
% of Electrical and Computer Engineering, Graduate School of Engineering,
% Air Force Institute of Technology, Wright Patterson Air Force Base, 1999.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Distance = min(pdist2(PopObj,PF),[],2);
    Score    = norm(Distance) / length(Distance);
end