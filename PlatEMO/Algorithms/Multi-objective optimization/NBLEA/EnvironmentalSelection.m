function Population = EnvironmentalSelection(Problem,Population,Offspring)
% The environmental selection of NBLEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    selected = randperm(length(Population),2);
    Pool     = [Population(selected),Offspring];
    [~,rank] = sort(CalFitness(Problem.C,Pool));
    Population(selected) = Pool(rank(1:length(selected)));
end