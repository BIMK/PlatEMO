function [Population, FrontNo, CrowdDis] = EnvironmentalSelection(Population, N)
% The environmental selection of PRDH
% Note: We assume the first objective is the selected feature ratio, and the
%       second objective is the classification error rate. 

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Ruwang Jiao

    % Remove duplicated solutions (search space)
    [~, index] = unique(Population.decs, 'rows');
    Population = Population(index);

    FrontNo = NDSort(Population.objs, inf);
    NextNo  = false(1, size(Population.objs, 1));
    
    % Directly save solutions in the first front
    NextNo(find(FrontNo==1)) = 1;
    NDpoints = Population(find(FrontNo==1));

    %% Duplication handling (objective space)
    for i = 2:max(FrontNo)
        FrontPop = [];
        No       = find(FrontNo==i);
        PopObj   = Population(No).objs;
        [~, ~, c] = unique(PopObj, 'rows', 'stable');
        for j=1:max(c)
            index = find(c==j);
            if size(index, 1) > 1
                selectedNo = DuplicationSelection(No(index), Population, NDpoints);
                FrontPop = [FrontPop, selectedNo];
            else
                FrontPop = [FrontPop, No(index)];
            end
        end
        NextNo(FrontPop) = 1;
    end
    Population = Population(NextNo);
    FrontNo    = FrontNo(NextNo);
    
    %% Problem reformulation and constraint handling
    % Find the solution with the minimum number of selected features
    index = findminobj1(Population.objs);
    PopObj = Population.objs;
    minobj2 = PopObj(index, 2);

    if sum(PopObj(:, 2)<=minobj2) >= N
        %% Selection among feasible solutions
        boolFeasible = PopObj(:, 2)<=minobj2;
        Population = Population(1:end, boolFeasible);
        FrontNo    = FrontNo(boolFeasible);
        sumt = 0;
        for MaxFNo=1:max(FrontNo)
            sumt = sumt + sum(FrontNo==MaxFNo);
            if sumt >= N
                break;
            end
        end
        Next = FrontNo < MaxFNo;
    
        %% Calculate the crowding distance of each solution
         CrowdDis = CrowdingDistance(Population.objs, FrontNo);
    
        %% Select the solutions in the last front based on their crowding distances
        Last     = find(FrontNo==MaxFNo);
        [~,Rank] = sort(CrowdDis(Last), 'descend');
        Next(Last(Rank(1:N-sum(Next)))) = true;
    
        %% Population for next generation
        Population = Population(Next);
        FrontNo    = FrontNo(Next);
        CrowdDis   = CrowdDis(Next);
    else
        %% Selection including infeasible solutions
        [~, rank]  = sort(PopObj(:, 1));
        if size(Population, 2) < N
            N = size(Population, 2);
        end
        Population = Population(rank(1:N));
        FrontNo = 1:N;
        CrowdDis = zeros(1, N);
    end
end

function selectedNo = DuplicationSelection(index, Population, NDpoints)
    %% Choose promising duplicated solutions (objective space) to survive
    NDobj = NDpoints.objs;
    NDdec = NDpoints.decs;
    Obj = Population(index).objs;
    Dec = Population(index).decs;
    c1 = find(min(abs(Obj(1,1)-NDobj(:,1)))==abs(Obj(1,1)-NDobj(:,1)));
    c2 = find(min(abs(Obj(1,2)-NDobj(:,2)))==abs(Obj(1,2)-NDobj(:,2)));
    o1 = [];
    for i=1:size(c1,1)
       tmp = pdist2(Dec(:,:), NDdec(c1(i,:), :), "hamming");
       o1 = [o1, tmp];
    end
    o1 = mean(o1, 2);
    o2 = [];
    for i=1:size(c2,1)
       index2 = find(NDdec(c2(i,:),:)==1);
       tmp = pdist2(Dec(:,index2), NDdec(c2(i,:), index2), "hamming");
       o2 = [o2, tmp];
    end
    o2 = mean(o2, 2);
    o  = [o1,o2];
    [FrontNO, ~] = NDSort(-o, inf);  
    No = find(FrontNO==1);
    selectedNo = index(No);
end

function minindex = findminobj1(obj)
    minindex = 1;
    for i=1:size(obj, 1)
        if obj(i, 1) < obj(minindex, 1) || (obj(i, 1)==obj(minindex, 1)&obj(i, 2)<obj(minindex, 2))
            minindex = i;
        end
    end
end