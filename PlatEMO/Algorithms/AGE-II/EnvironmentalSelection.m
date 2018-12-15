function Population = EnvironmentalSelection(Population,ArcObj,N)
% The environmental selection of AGE-II

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    NP = length(Population);
    NA = size(ArcObj,1);
    
    %% Discard the offsprings dominated by the increment of archive
    discard = false(1,NP);
    for i = N+1 : NP
        discard(i) = any(all(repmat(Population(i).obj,NA,1)>=ArcObj+1,2));
    end
    Population(discard) = [];
    NP = length(Population);
    
    %% Remove solutions by fast approximation
    alpha = zeros(NP,NA);
    for i = 1 : NA
        % The formula in the original paper is incorrect, which has been
        % revised here
        alpha(:,i) = max(Population.objs-repmat(ArcObj(i,:),NP,1),[],2);
    end
    [rho,rank] = sort(alpha,1);
    % Delete the solution one by one
    Remain = 1 : NP;
    while length(Remain) > N
        % Calculate the approximations when each solution is eliminated
        % from the population
        S = zeros(length(Remain),NA);
        for i = 1 : length(Remain)
            temp = rank(1,:) == Remain(i);
            S(i,~temp) = rho(1,~temp);
            S(i,temp)  = rho(2,temp);
        end
        % Delete the worst solution in the population and update the
        % variables
        [~,worst] = sortrows(sort(S,2,'descend'));
        remain = rank ~= Remain(worst(1));
        rho    = reshape(rho(remain),length(Remain)-1,NA);
        rank   = reshape(rank(remain),length(Remain)-1,NA);
        Remain(worst(1)) = [];
    end
    % Population for next generation
    Population = Population(Remain);
end