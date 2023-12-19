function [NearStage,Fitness,FitnessLayer,LayerMax] = ControlStage(SparseRate,NearStage,Mask,Dec,Fitness,FitnessLayer,LayerMax,Problem)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Stage = ceil(Problem.FE/(Problem.maxFE/10));
    if (Stage ~= NearStage)
        NearStage = Stage;
        [FitnessLayer,LayerMax] = UpdateLayer(SparseRate,Stage,Fitness,Problem,Mask);
    end        
end