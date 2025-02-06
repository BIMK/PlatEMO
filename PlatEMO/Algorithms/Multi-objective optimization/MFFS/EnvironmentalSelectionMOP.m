function [Population, FrontNo, CrowdDis] = EnvironmentalSelectionMOP(Population, N)
% Environmental selection of MFFS for the multi-objective task

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Ruwang Jiao

    % Calculate average distance of the parent population in the search space
    PopDec = Population(1:N).decs;
    aDist  = pdist2(round(PopDec), round(PopDec), "hamming");
    aveDis = sum(sum(tril(aDist, 0)))./sum(1:(size(PopDec,1)-1));
    
    % Remove duplicated solutions in the population
    PopDec     = Population.decs;
    [~, index] = unique(PopDec, 'rows');
    Population = Population(index);
    PopObj     = Population.objs;
    PopDec     = Population.decs;
    
    % Nondominated sorting
    FrontNo = NDSort(PopObj, inf);
    NextNo  = false(1, size(PopObj, 1));
    
    % Directly save solutions in the first front
    First     = find(FrontNo==1);
    NextNo(First) = 1;
    NDpoints  = Population(First);
    
    while sum(NextNo) < N || sum(NextNo)==size(PopObj, 1)
        for i = 2:max(FrontNo)
            FrontPop = [];
            No       = find(FrontNo==i);
            No       = No(NextNo(No)==0);
            Pop      = Population(No);
            PopObj   = Pop.objs;
            [Obj, ~, c] = unique(PopObj, 'rows', 'stable');
            for j=1:max(c)
                index = find(c==j);
                if size(index, 1) > 1
                    selectedNo = DuplicationSelection(No(index), Population, NDpoints, NextNo, aveDis);
                    FrontPop = [FrontPop, selectedNo];
                else
                    FrontPop = [FrontPop, No(index)];
                end
            end
            NextNo(FrontPop) = 1;
            if sum(NextNo) >= N
                break;
            end
        end
    end
    
    Population = Population(NextNo);
    FrontNo    = FrontNo(NextNo);
    MaxFNo     = max(FrontNo);
    Next       = FrontNo < MaxFNo;
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(Population.objs, FrontNo);
    
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo == MaxFNo);
    [~,Rank] = sort(CrowdDis(Last), 'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
end
    
 function selectedNo = DuplicationSelection(index, Population, NDpoints, Next, aveDis)
     % Choose one solution from many duplicated solutions (objective space)
     NDobj = NDpoints.objs;
     NDdec = NDpoints.decs;
     Nextobj = Population(Next).objs;
     Nextdec = Population(Next).decs;
     Obj = Population(index).objs;
     Dec = Population(index).decs;
     [r, c] = ismember(Obj(1,1), Nextobj(:,1));
     if r == 1
         aDist = pdist2(Dec, Nextdec(c,:), "hamming");
         if sum(index(aDist>=aveDis)) >= 1
             selectedNo = index(aDist>=aveDis);
         else
            [~, No] = max(sum(aDist,2), [], 1);
            selectedNo = index(No);
         end
     else
         Dis = abs(repmat(Obj(1,1),size(NDobj,1),1) - NDobj(:,1));
         [~, no] = min(Dis, [], 1);
         aDist = pdist2(Dec, NDdec(no,:), "hamming");
         if sum(index(aDist>=aveDis)) >= 1
              selectedNo = index(aDist>=aveDis);
         else
            [~, No] = max(sum(aDist,2), [], 1);
            selectedNo = index(No);
         end
     end
 end