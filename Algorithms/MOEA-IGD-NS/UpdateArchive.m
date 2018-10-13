function Archive = UpdateArchive(Population,NA)
% Update the archive in MOEA/IGD-NS

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    %% Detect the non-dominated solutions
    Population = Population(NDSort(Population.objs,1)==1);
    
    %% Select the extreme solutions
    Choose          = false(1,length(Population)); 
    [~,extreme]     = max(Population.objs,[],1);
    Choose(extreme) = true;
    
    %% Select other solutions by truncation
    if sum(Choose) > NA
        selected = find(Choose);
        Choose   = selected(randperm(length(selected),NA));
    else
        Cosine = 1 - pdist2(Population.objs,Population.objs,'cosine');
        Cosine(logical(eye(length(Cosine)))) = 0;
        while sum(Choose) < NA && ~all(Choose)
            unSelected = find(~Choose);
            [~,x]      = min(max(Cosine(~Choose,Choose),[],2));
            Choose(unSelected(x)) = true;
        end
    end
    Archive = Population(Choose);
end