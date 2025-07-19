function newPopulation = addNoise(init_population, Nini, n)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Num           = size(init_population,1);
    newPopNumber  = Nini - Num;
    newPopulation = zeros(newPopNumber,n);  
    for i = 1 : newPopNumber
        index  = randperm(Num,1);
        newPopulation(i,:) = init_population(index,:);
        index2 = randperm(n,round(n*0.3));
        temp   = init_population(index,index2) + normrnd(0,0.5,[1 round(n*0.3)]);
        upp    = temp>1;
        loww   = temp<0;
        temp(upp)  = (init_population(index,upp)+1)/2;
        temp(loww) = (init_population(index,loww)-0)/2;
        newPopulation(i,index2) = temp;
    end
end