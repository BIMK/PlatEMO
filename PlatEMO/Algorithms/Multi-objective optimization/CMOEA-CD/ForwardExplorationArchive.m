function FA = ForwardExplorationArchive(FA, Offspring, zmin, Ns, Add)
% The Forward Exploration Archive of CMOEA-CD
% Add = 1 --- environmental selection of SPEA2
% Add = 2 --- environmental selection of NSGA-II
% Add = 3 --- environmental selection of modified NSGA-III

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhe Liu

    %% Dominance relation calculation
    FA = [FA, Offspring];
    NonDominated = DominationCal(FA, 0);
    FA = FA(NonDominated);
    PopObj = FA.objs;
    [P, M] = size(PopObj);  
    if P <= Ns
        FA = FA([1:P, unidrnd(P,[1, Ns - P])]);
    elseif Add == 1
        zmax   = max(PopObj, [], 1);        
        PopObj = (PopObj - zmin) ./ (zmax - zmin);
        Del    = Truncation(PopObj, P - Ns);
        remained_solution = 1: P;
        remained_solution(Del) = [];
        FA = FA(remained_solution);
    elseif Add == 2
        zmax     = max(PopObj, [], 1);        
        PopObj   = (PopObj - zmin) ./ (zmax - zmin);
        CrowdDis = CrowdingDistance(PopObj);
        [~,Rank] = sort(CrowdDis,'descend');  
        FA = FA(Rank(1:Ns));       
    elseif Add == 3     
        [W, Ns]  = UniformPoint(Ns, M);     
        zmax     = max(PopObj, [], 1); 
        PopObj   = (PopObj - zmin) ./ (zmax - zmin);
        Distance = sqrt(sum(PopObj.^2, 2));
        Angle_Pop_to_W = sin(acos(1 - pdist2(W, PopObj, 'cosine')));
        S  = FA;
        FA = [];
        for i = 1 : Ns
            Fitness   = Distance' .* Angle_Pop_to_W(i,:);
            [~,index] = min(Fitness);
            FA = [FA,S(index)];
        end
    end   
end

function Del = Truncation(PopObj,K)
% Select part of the solutions by truncation

    %% Truncation
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,size(PopObj,1));
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end