function [Population,FrontNo] = Reduce(Population,FrontNo)
% Delete one solution from the population

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Identify the solutions in the last front
    FrontNo   = UpdateFront(Population.objs,FrontNo);
    LastFront = find(FrontNo==max(FrontNo));
    PopObj    = Population(LastFront).objs;
    [N,M]     = size(PopObj);
    
    %% Calculate the contribution of hypervolume of each solution
    deltaS = inf(1,N);
    if M == 2
        [~,rank] = sortrows(PopObj);
        for i = 2 : N-1
            deltaS(rank(i)) = (PopObj(rank(i+1),1)-PopObj(rank(i),1)).*(PopObj(rank(i-1),2)-PopObj(rank(i),2));
        end
    elseif N > 1
        deltaS = CalHV(PopObj,max(PopObj,[],1)*1.1,1,10000);
    end
    
    %% Delete the worst solution from the last front
    [~,worst] = min(deltaS);
    FrontNo   = UpdateFront(Population.objs,FrontNo,LastFront(worst));
    Population(LastFront(worst)) = [];
end