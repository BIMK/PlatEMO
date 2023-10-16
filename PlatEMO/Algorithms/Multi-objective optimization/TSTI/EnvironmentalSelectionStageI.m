function [Population,fitness] = EnvironmentalSelectionStageI(Population,PopObj,OffObj,N)
% The environmental selection of TiGE_1

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

	[FrontNo,~] = NDSort([PopObj;OffObj],N);
    fcv         = Calculate_fcv(Population);
    fitness     = FrontNo' + fcv./(fcv+1);
    [~,index]   = sort(fitness);
    Population  = Population(index(1:N));
    fitness     = fitness(index(1:N));
end