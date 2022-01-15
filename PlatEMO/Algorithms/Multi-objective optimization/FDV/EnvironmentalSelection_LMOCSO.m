function Population = EnvironmentalSelection_LMOCSO(Population,V,theta)
% The environmental selection of LMOCSO

%  Copyright (C) 2021 Xu Yang
%  Xu Yang <xuyang.busyxu@qq.com> or <xuyang369369@gmail.com>

    Population = Population(NDSort(Population.objs,1)==1);%取非支配排序第一层
    PopObj = Population.objs;
    [N,M]  = size(PopObj);
    NV     = size(V,1);
    
    %% Translate the population
    PopObj = PopObj - repmat(min(PopObj,[],1),N,1);
    
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
    Population = Population(Next(Next~=0));
end