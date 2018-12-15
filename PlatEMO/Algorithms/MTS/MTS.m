function MTS(Global)
% <algorithm> <M>
% Multiple trajectory search
% popsize           --- 40 --- Size of the population
% ofLocalSearchTest ---  5 --- Number of iterations for determining the best local search
% ofLocalSearch     --- 45 --- Number of iterations for applying the best local search
% ofForeground      ---  5 --- Number of best solutions for local search

%------------------------------- Reference --------------------------------
% L. Y. Tseng and C. Chen, Multiple trajectory search for unconstrained /
% constrained multi-objective optimization, Proceedings of the IEEE
% Congress on Evolutionary Computation, 2009, 1951-1958.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [popsize,ofLocalSearchTest,ofLocalSearch,ofForeground] = Global.ParameterSet(40,5,45,5);

    %% Generate random population
    % Generate the SOA
    [~,SOA] = sort(rand(popsize,Global.D));
    % Map the SOA to the decision space
    Decs = SOA./popsize.*repmat(Global.upper-Global.lower,popsize,1) + repmat(Global.lower,popsize,1);
    % Generate the initial solutions
    Population  = INDIVIDUAL(Decs);
    AppSet      = [];
    Enable      = true(1,popsize);
    Improve     = true(1,popsize);
    Grade       = zeros(1,popsize);
    SearchRange = repmat((Global.upper-Global.lower)/2,Global.N,1);

    %% Optimization
    while Global.NotTermination(AppSet)
        for i = find(Enable)
            Grade(i)     = 0;
            LS_TestGrade = zeros(1,3);
            for j = 1 : ofLocalSearchTest
                [grade1,Population(i),SearchRange(i,:),Improve(i),AppSet] = LocalSearch1(Global,Population(i),SearchRange(i,:),Improve(i),AppSet);
                [grade2,Population(i),SearchRange(i,:),Improve(i),AppSet] = LocalSearch2(Global,Population(i),SearchRange(i,:),Improve(i),AppSet);
                [grade3,Population(i),SearchRange(i,:),Improve(i),AppSet] = LocalSearch3(Global,Population(i),SearchRange(i,:),Improve(i),AppSet);
                LS_TestGrade(1) = LS_TestGrade(1) + grade1;
                LS_TestGrade(2) = LS_TestGrade(2) + grade2;
                LS_TestGrade(3) = LS_TestGrade(3) + grade3;
            end
            [~,best]    = max(LS_TestGrade);
            LocalSearch = {@LocalSearch1,@LocalSearch2,@LocalSearch3};
            for j = 1 : ofLocalSearch
                [grade,Population(i),SearchRange(i,:),Improve(i),AppSet] = LocalSearch{best}(Global,Population(i),SearchRange(i,:),Improve(i),AppSet);
                Grade(i) = Grade(i) + grade;
            end
        end
        [~,rank]  = sort(Grade,'descend');
        Enable(:) = false;
        Enable(rank(1:ofForeground)) = true;
        if Global.evaluated >= Global.evaluation
            AppSet = AdjustAppSet(AppSet,Global.N);
        end
    end
end