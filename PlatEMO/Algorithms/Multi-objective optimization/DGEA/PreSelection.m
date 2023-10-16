function [Population, FrontNo] = PreSelection(Population,V,theta,RefNo)
% The preselection strategy in DGEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He

    %% Preseleting some well-converged solutions
    NV     = size(V,1);
    [Front,MaxFront] = NDSort(Population.objs,min([NV,length(Population)]));
    if  sum(Front==1) >= RefNo  
        Next0 = Front < MaxFront;
        Last = find(Front == MaxFront);
    else
        Next0 = Front == 1;
        Last = find(Front >1);
    end

    %% Selecting those diversity-related but well-converged solutions
    PopObj = Population(Last).objs;
    [N,M]  = size(PopObj);

    %% Translate the population
    PopObj = PopObj - repmat(min(Population.objs,[],1),N,1);
    
    %% Calculate the smallest angle value between each vector and others
    cosine = 1 - pdist2(V,V,'cosine');
    cosine(logical(eye(length(cosine)))) = 0;
    gamma  = min(acos(cosine),[],2);

    %% Associate each solution to a reference vector
    Angle = acos(1-pdist2(PopObj,V,'cosine'));
    [~,associate] = min(Angle,[],2);

    %% Select one solution for each reference vector
    Next = zeros(1,NV);
    for i = unique(associate)'
        current = find(associate==i);
        % Calculate the APD value of each solution
        APD = (1+M*theta*Angle(current,i)/gamma(i)).*sqrt(sum(PopObj(current,:).^2,2));
        % Select the one with the minimum APD value
        [~,best] = min(APD);
        Next(i) = current(best);
    end
    % Population for next generation
    Next0(Last(Next(Next~=0))) = true;
    Population = Population(Next0);
    FrontNo = Front(Next0);
end