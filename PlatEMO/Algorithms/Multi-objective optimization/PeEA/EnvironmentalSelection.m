function Population = EnvironmentalSelection(Population,T)
% The environmental selection of PeEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Li Li

    [Fitness,extreme] = CalFitness(Population.objs);
    extreme = unique(extreme); %eliminate the identical solutions
    [~,mm] = size(extreme);  
   
	%% Calculate the angle between each two solutions
    Angle = acos(1-pdist2(Population.objs,Population.objs,'cosine'));
    Angle(logical(eye(length(Population)))) = inf;

    %% Angle-based tournament selection
    Remain = 1 : length(Population);
    Remain(extreme) = []; %eliminate extremes first
    while length(Remain) > (T-mm) %mins mm first
        % Identify the two solutions A and B with the minimum angle
        [sortA,rank1] = sort(Angle(Remain,Remain),2); %sortA is the sorted angle, rank1 NXN index of matric angle
        [~,rank2]     = sortrows(sortA);              %rank2 Nx1 index of sorted angle
        A = rank2(1);                                 %A, with minimum angle
        B = rank1(A,1);                               %B corrponding to A with minimum angle
        % Eliminate one of A and B
        if  Fitness(Remain(A)) > Fitness(Remain(B))
            Remain(A) = [];
        else
            Remain(B) = [];
        end
    end
    % Population for next generation
    Population = Population([Remain,extreme]); %add extremes into pop
end