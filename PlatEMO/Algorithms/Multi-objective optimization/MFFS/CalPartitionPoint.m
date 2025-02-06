function PartitionSet = CalPartitionPoint(Population, FrontNo, Pset)
% Calculate partion point between two largest nondominated soltuions
    
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Ruwang Jiao

    [PopObj, ic] = unique(Population.objs, 'rows');
    FrontNo      = FrontNo(ic);
    CrowdAngle   = zeros(size(PopObj, 1), size(PopObj, 1));
    Fronts       = setdiff(unique(FrontNo), inf);
    Front        = find(FrontNo==Fronts(1));
    Fmin         = min(PopObj(Front, :), [], 1);
    [~, Rank] = sortrows(PopObj(Front, 1));
    for j = 1 : length(Front) - 1
        CrowdAngle(Front(Rank(j)), Front(Rank(j+1))) = acos(1-pdist2((PopObj(Front(Rank(j+1)),:) - Fmin), (PopObj(Front(Rank(j)),:) - Fmin), 'cosine'));
    end
    
    if length(Front)==1
        PartitionSet = zeros(2, 2);
    elseif length(Front)==2
        [r, c] = find(CrowdAngle == max(max(CrowdAngle)));
        PartitionSet = [(PopObj(r,:) + PopObj(c,:))./2; (PopObj(r,:) + PopObj(c,:))./2];
        CrowdAngle(r, c) = 0;
    else
        [r, c] = find(CrowdAngle == max(max(CrowdAngle)));
        MidP1 = (PopObj(r,:) + PopObj(c,:))./2;
        CrowdAngle(r, c) = 0;
        [r, c] = find(CrowdAngle == max(max(CrowdAngle)));
        MidP2 = (PopObj(r,:) + PopObj(c,:))./2;
        CrowdAngle(r, c) = 0;
        PartitionSet = [MidP1; MidP2];
    end
    if ismember(PartitionSet(1,:), Pset, 'rows') 
        if max(max(CrowdAngle))~=0
            [r, c] = find(CrowdAngle == max(max(CrowdAngle)));
            PartitionSet(1,:) = (PopObj(r,:) + PopObj(c,:))./2;
            CrowdAngle(r, c) = 0;
        else
            PartitionSet(1,:) = PartitionSet(2,:);
        end
    end
    if ismember(PartitionSet(2,:), Pset, 'rows')
        if max(max(CrowdAngle))~=0
            [r, c] = find(CrowdAngle == max(max(CrowdAngle)));
            PartitionSet(2,:) = (PopObj(r,:) + PopObj(c,:))./2;
        else
            PartitionSet(2,:) = PartitionSet(1,:);
        end
    end
end