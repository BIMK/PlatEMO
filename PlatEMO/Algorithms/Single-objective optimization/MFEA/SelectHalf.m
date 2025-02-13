function  Parent = SelectHalf(SubPopulation)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Parent = {};
    for i = 1 : size(SubPopulation,2)
        pop       = SubPopulation{i};
        [~,rank]  = sort(pop.objs);
        Parent{i} = pop(rank(1:floor(size(pop,2)/2)));
    end
end