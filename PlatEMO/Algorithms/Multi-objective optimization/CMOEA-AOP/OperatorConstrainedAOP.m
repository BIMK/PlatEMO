function Offspring = OperatorConstrainedAOP(Problem, Population, MatingPool, action, i)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    rates    = action;
    ratesNum = floor(rates.*Problem.N);
    for item = 1 : length(action)
        if ratesNum(item) < 4
            ratesNum(item) = 4; % Set the minimum value to 4 so that each operator can be called normally
        end
    end
    MatingPool = repmat(MatingPool,1,20);
    Offspring1 = OperatorGAhalf(Problem,Population{i}(MatingPool(1:ratesNum(1)*2)));
    Offspring2 = OperatorDE(Problem,Population{i}(randi(ceil(Problem.N),1,ratesNum(2))),Population{i}(randi(ceil(Problem.N),1,ratesNum(2))),Population{i}(randi(ceil(Problem.N),1,ratesNum(2))));
    p3         = Population{i};
    Offspring3 = DEBest(Problem, p3(randi(ceil(Problem.N),1,ratesNum(3))),ratesNum(3));
    Offspring  = [Offspring1 Offspring2 Offspring3];
end