function [SparseRate,Fitness] = MGCEA_ClusterLabel(FitnessInit,Fitness)
    Num1 = sum(Fitness == 1);
    Num2 = sum(Fitness == 2);
    Num1Value = sum(FitnessInit(Fitness == 1));
    Num2Value = sum(FitnessInit(Fitness == 2));
    if Num1Value < Num2Value
        SparseRate = Num1/(Num1 + Num2);
        Fitness(Fitness == 1) = 11;
        Fitness(Fitness == 2) = 12;
    else   
        SparseRate = Num2/(Num1 + Num2);
        Fitness(Fitness == 1) = 12;
        Fitness(Fitness == 2) = 11;
    end
end