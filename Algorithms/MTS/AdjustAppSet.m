function AppSet = AdjustAppSet(AppSet,N)
% Reduce the size of approximation set

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    AppObj = AppSet.objs;
    
    %% Selete the solutions having the smallest values on at least one objective
    Choose = false(1,length(AppSet));
    for i = 1 : size(AppObj,2)
        Choose(AppObj(:,i)==min(AppObj(:,i))) = true;
    end
    
    %% Select other solutions one by one
    if sum(Choose) > N
        selected = find(Choose);
        Choose   = selected(randperm(length(selected),N));
    else
        Distance = pdist2(AppObj,AppObj);
        Distance(logical(eye(length(Distance)))) = inf;
        while sum(Choose) < N && ~all(Choose)
            unSelected = find(~Choose);
            [~,x]      = max(min(Distance(~Choose,Choose),[],2));
            Choose(unSelected(x)) = true;
        end
    end
    AppSet = AppSet(Choose);
end