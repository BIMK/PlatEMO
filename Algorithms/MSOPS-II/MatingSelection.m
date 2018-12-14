function Parents = MatingSelection(Population,Archive)
% The mating selection of MSOPS-II

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N  = length(Population);
    NA = length(Archive);

    %% Select one solution in archive for each solution in population
    Parent1 = randi(NA,1,N);
    Parent2 = randi(NA,1,N);
    % Use binary tournament selection to select one solution in archive
    % which is nearer to the solution in population
    temp    = rand(1,N) < 0.5;
    Choose1 = true(N,1);
    Choose1(temp)  = sum((Archive(Parent1(temp)).objs-Population(temp).objs).^2,2) < sum((Archive(Parent2(temp)).objs-Population(temp).objs).^2,2);
    Choose1(~temp) = sum((Archive(Parent1(~temp)).decs-Population(~temp).decs).^2,2) < sum((Archive(Parent2(~temp)).decs-Population(~temp).decs).^2,2);
    MatingPool = [Parent1(Choose1),Parent2(~Choose1)];
    
    %% Use each solution in population and each solution in mating pool to generate one offspring
    Parents = [Population,Archive(MatingPool)];
end