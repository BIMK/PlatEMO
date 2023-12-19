function [FitnessLayer,LayerMax] = UpdateLayer(SparseRate,Stage,Fitness,Problem,Mask)
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