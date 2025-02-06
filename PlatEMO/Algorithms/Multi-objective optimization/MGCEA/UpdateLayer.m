function [FitnessLayer,LayerMax] = UpdateLayer(SparseRate,Stage,Fitness,Problem,Mask)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    GroupNum = ceil(11 - Stage)/100*Problem.D;
    GroupNum = ceil(SparseRate*10*1*GroupNum);
    if sum(sum(Mask)) == 0
        [~,FitnessIndex] = sort(Fitness + rand(1,Problem.D));
    else
        [~,FitnessIndex] = sort(Fitness + sum(Mask == 0)./100000);
    end
    FitnessIndexLayer = ceil((1:Problem.D)./GroupNum);
    FitnessLayer = zeros(1,Problem.D);
    FitnessLayer(FitnessIndex) = FitnessIndexLayer;
    LayerMax = max(FitnessLayer);
end