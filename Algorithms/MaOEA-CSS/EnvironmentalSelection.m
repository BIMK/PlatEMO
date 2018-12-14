function Population = EnvironmentalSelection(Population,Zmin,t,K)
% The environmental selection of MaOEA-CSS

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Calculate the distance between each solution to the ideal point
    PopObj = Population.objs - repmat(Zmin,length(Population),1);
    Con    = sqrt(sum(PopObj.^2,2));

	%% Calculate the angle between each two solutions
    Angle = acos(1-pdist2(PopObj,PopObj,'cosine'));
    Angle(logical(eye(length(Population)))) = inf;
    
    %% Eliminate solutions one by one
    Remain = 1 : length(Population);
    while length(Remain) > K
        % Identify the two solutions A and B with the minimum angle
        [sortA,rank1] = sort(Angle(Remain,Remain),2);
        [~,rank2]     = sortrows(sortA);
        A = rank2(1);
        B = rank1(A,1);
        % Eliminate one of A and B
        if Con(Remain(A)) - Con(Remain(B)) > t
            Remain(A) = [];
        elseif Con(Remain(B)) - Con(Remain(A)) > t
            Remain(B) = [];
        else
            Remain(A) = [];
        end
    end
    % Population for next generation
    Population = Population(Remain);
end