function [SubPopulation,Rank] = Sort(Problem,SubPopulation)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Rank = {};
    for i = 1 : length(Problem.SubD)
        [~,FrontNo,CrowdDis] = EnviSelect(SubPopulation{i},length(SubPopulation{i}));
        [~,rank]             = sortrows([FrontNo',-CrowdDis']);
        SubPopulation{i}     = SubPopulation{i}(rank);
        Rank{i}              = 1 : length(rank);
    end
end