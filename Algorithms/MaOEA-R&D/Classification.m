function [Subpopulation,Z] = Classification(Population,W)
% Classify the population into each subpopulation

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Calculate ASF between each solution and classification vector
    PopObj = Population.objs;
    [N,M]  = size(PopObj);
    Z      = min(PopObj,[],1);
    ASF    = zeros(N,M);
    for i = 1 : M
        ASF(:,i) = max((PopObj-repmat(Z,N,1))./repmat(W(i,:),N,1),[],2);
    end
    
    %% Classification
    class    = zeros(N,1);
    external = true(1,N);
    lacking  = true(1,M);
    while any(lacking)
        % Redistribute the external solutions
        lacks = find(lacking);
        [~,subclass]    = min(ASF(external,lacks),[],2);
        class(external) = lacks(subclass);
        external = false(1,N);
        % Clip the exceeded subpopulations
        for i = 1 : M
            sub = find(class==i);
            if length(sub) > N/M
                [~,rank] = sort(ASF(sub,i));
                external(sub(rank(N/M+1:end))) = true;
            end
            lacking(i) = length(sub) < N/M;
        end
    end
    Subpopulation = cell(1,M);
    for i = 1 : M
        Subpopulation{i} = Population(class==i);
    end
end