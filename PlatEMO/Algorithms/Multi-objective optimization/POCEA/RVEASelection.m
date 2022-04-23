function Population = RVEASelection(Population,V,popsize,theta)
% Uniformity optimization by RVEA

% Copyright (c) 2020-2021 Cheng He

    PopObj = Population.objs;
    [N,~]  = size(PopObj);
    NV     = size(V,1);
    
    %% Translate the population
    PopObj = PopObj - repmat(min(PopObj,[],1),N,1);
    
    %% Calculate the degree of violation of each solution
    CV = sum(max(0,Population.cons),2);
    
    %% Calculate the smallest angle value between each vector and others
    cosine = 1 - pdist2(V,V,'cosine');
    cosine(logical(eye(length(cosine)))) = 0;
    gamma  = min(acos(cosine),[],2);

    %% Associate each solution to a reference vector
    Angle = acos(1-pdist2(PopObj,V,'cosine'));
    [~,associate] = min(Angle,[],2);

    %% Select one solution for each reference vector
    Next = zeros(1,NV);
    pf = sum(CV(randi(N,[N,1]))<1e-6)/N;
    for i = unique(associate)'
        cv = CV(associate==i);
        Ns = sum(associate==i);
        subN = ceil(popsize/numel(unique(associate)'));
		% Ensure the subpopulation size
        if Ns < subN
            epsilon = max(cv);
        else
            epsilon = min(cv)*(1-pf)+mean(cv)*pf;
        end
        current1 = find(associate==i & CV<=epsilon);
        current2 = find(associate==i & CV>epsilon);
        if ~isempty(current1)
            % Calculate the APD value of each solution
            APD = (1+theta*Angle(current1,i)/gamma(i)).*sqrt(sum(PopObj(current1,:).^2,2));
            [~,best] = min(APD);
            Next(i)  = current1(best);
        elseif ~isempty(current2)
            % Select the one with the minimum CV value
            [~,best] = min(CV(current2));
            Next(i)  = current2(best);
        end
    end
    % Population for next generation
    Population = Population(Next(Next~=0));
end