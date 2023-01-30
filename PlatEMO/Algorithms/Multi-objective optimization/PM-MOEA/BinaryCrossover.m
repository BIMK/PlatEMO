function Offspring = BinaryCrossover(Parent1,Parent2)
% Unbalanced binary crossover

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Offspring = logical(Parent1);
    for i = 1 : size(Offspring,1)
        diff = find(Parent1(i,:)~=Parent2(i,:));
        MOne = mean(Offspring(i,diff));
        r    = min(min(0.5,2*MOne),2*(1-MOne));
        rate = zeros(1,length(diff));
        rate(Offspring(i,diff))  = r/2/MOne;
        rate(~Offspring(i,diff)) = r/2/(1-MOne);
        exchange                    = rand(1,length(diff)) < rate;
        Offspring(i,diff(exchange)) = ~Offspring(i,diff(exchange));
    end
end