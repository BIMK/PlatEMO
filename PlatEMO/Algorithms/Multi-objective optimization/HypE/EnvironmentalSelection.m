function Population = EnvironmentalSelection(Population,N,RefPoint,nSample)
% The environmental selection of HypE

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

    %% Select the solutions in the last front
    Last   = find(FrontNo==MaxFNo);
    Choose = true(1,length(Last));
    while sum(Choose) > N-sum(Next)
        drawnow('limitrate');
        Remain  = find(Choose);
        F       = CalHV(Population(Last(Remain)).objs,RefPoint,sum(Choose)-N+sum(Next),nSample);
        [~,del] = min(F);
        Choose(Remain(del)) = false;
    end
    Next(Last(Choose)) = true;
    % Population for next generation
    Population = Population(Next);
end