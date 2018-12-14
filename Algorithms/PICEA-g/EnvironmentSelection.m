function [Population,Goal] = EnvironmentSelection(Population,Goal,N)
% The environmental selection of PICEA-g

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    NP    = length(Population);
    NGoal = size(Goal,1);

    %% Calculate ng
    FdG = false(NP,NGoal);
    ng  = zeros(1,NGoal);
    for i = 1 : NP
        x = all(repmat(Population(i).obj,NGoal,1)-Goal<=0,2);
        FdG(i,x) = true;
        ng(x)    = ng(x)+1;
    end
    
    %% Calculate Fs
    Fs = zeros(1,NP);
    for i = 1 : NP
        Fs(i) = sum(1./ng(FdG(i,:)));
    end
    
    %% Calculate Fg
    Fg = zeros(1,NGoal);
    for i = 1 : NGoal
        if ng(i) == 0
            Fg(i) = 0.5;
        else
            Fg(i) = 1/(1+(ng(i)-1)/(NP-1));
        end
    end   
    
    %% Select half of the solutions
    ND = find(NDSort(Population.objs,1)==1);
    if length(ND) < N
        Fs(ND) = inf;
        [~,Rank] = sort(Fs,'descend');
        Next = Rank(1:N);
    else
        [~,Rank] = sort(Fs(ND),'descend');
        Next = ND(Rank(1:N));
    end
    Population = Population(Next);
    
    %% Select half of the goals
    [~,Rank] = sort(Fg,'descend');
    Next = Rank(1:N);
    Goal = Goal(Next,:);
end