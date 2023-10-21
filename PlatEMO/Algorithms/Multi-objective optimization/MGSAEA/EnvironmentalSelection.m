function [PopDec,PopObj,Fitness] = EnvironmentalSelection(PopDec,PopObj,NI,M,status)
% Environmental selection

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Delete the duplicated points
    [~, Unduplicated] = unique(PopObj(:,1:M),'rows');
    PopDec            = PopDec(Unduplicated,:);
    PopObj            = PopObj(Unduplicated,:);
    
    %% Calculate the fitness of each solution
    if nargin == 4
        Fitness = CalFitness(PopObj);
    else
        if status == 1
            RealObj = PopObj(:,1:M);
            CV      = PopObj(:,end);           
            Fitness = CalFitness(RealObj,max(0,CV)); 
        elseif status == 2
            RealObj = PopObj(:,1:M);
            PopCon  = PopObj(:,M+1:end);
            CV      = sum(max(0,PopCon),2);            
            Fitness = CalFitness([RealObj,CV]);
        elseif status == 3
            Fitness = CalFitness(PopObj);
        end
    end
    
    %% Environmental selection
    if nargin == 4
        Next = Fitness < 1;
        if sum(Next) < NI
            [~,Rank] = sort(Fitness);
            Next(Rank(1:NI)) = true;
        elseif sum(Next) > NI
            Del  = Truncation(PopObj(Next,:),sum(Next)-NI);
            Temp = find(Next);
            Next(Temp(Del)) = false;
        end
    else
        Next = Fitness < 1;
        if sum(Next) < NI
            [~,Rank] = sort(Fitness);
            Next(Rank(1:NI)) = true;
        elseif sum(Next) > NI
            if status ~=3
                RealObj = PopObj(:,1:M);
            else
                RealObj = PopObj;
            end
            Del  = Truncation(RealObj(Next,:),sum(Next)-NI);
            Temp = find(Next);
            Next(Temp(Del)) = false;
        end
    end
    PopDec     = PopDec(Next,:);
    PopObj     = PopObj(Next,:);   
    Fitness    = Fitness(Next);
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