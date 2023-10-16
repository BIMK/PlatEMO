function Population = EnvironmentalSelection(Population,N)
% The environmental selection of RM-MEDA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,N);
    Next = FrontNo < MaxFNo;

    %% Delete the solutions in the last front by crowding distance
    Last = find(FrontNo==MaxFNo);
    while length(Last) > N - sum(Next)
        [~,worst]   = min(CrowdingDistance(Population(Last).objs));
        Last(worst) = [];
    end
    Next(Last) = true;
    % Population for next generation
    Population = Population(Next);
end