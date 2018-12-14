function Population = EnvironmentalSelection(Population,Zmin,N,K,L)
% The environmental selection of MaOEA-DDFC

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,N);
    Next = FrontNo < MaxFNo;

    %% Select the solutions in the last front
    Last   = find(FrontNo==MaxFNo);
    Choose = LastSelection(Population(Next).objs,Population(Last).objs,Zmin,N,K,L);
    Next(Last(Choose)) = true;
    % Population for next generation
    Population = Population(Next);
end

function Choose = LastSelection(P,F,Zmin,Total,K,L)
% Select part of the solutions in the last front

    PopObj = [P;F];
    [N,M]  = size(PopObj);

    %% Calculate the FC value
    FC = CalFC(PopObj,Zmin);

    %% Projection
    % Identify the ideal point
    Zmin = min(PopObj,[],1);
    % Identify the extreme points
    W = zeros(M) + 1e-6;
    W(logical(eye(M))) = 1;
    ASF = zeros(N,M);
    for i = 1 : M
        ASF(:,i) = max((PopObj-repmat(Zmin,N,1))./repmat(W(i,:),N,1),[],2);
    end
    [~,extreme] = min(ASF,[],1);
    % Calculate the intercepts
    Hyperplane = PopObj(extreme,:)\ones(M,1);
    a = (1./Hyperplane)';
    if any(isnan(a))
        a = max(PopObj,[],1);
    end
    % Normalization
    PopObj = (PopObj-repmat(Zmin,N,1))./repmat(a-Zmin,N,1);
    % Projection
    PopObj = PopObj./(repmat(sum(PopObj,2),1,M));
    
    %% Select the solutions in the last front one by one
    Choose = false(1,N);
    Choose(1:size(P,1)) = true;
    % The distance between each two solutions
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    while sum(Choose) < Total
        % Direction-based selection
        if ~any(Choose)
            Dis = sort(Distance(~Choose,~Choose),2);
        else
            Dis = sort(Distance(~Choose,Choose),2);
        end
        DD       = sum(1./(Dis(:,1:min(K,end))),2);
        [~,rank] = sort(DD);
        Remain   = find(~Choose);
        R        = Remain(rank(1:min(L,end)));
        % Convergence-based roulette-wheel selection
        Choose(R(RouletteWheelSelection(1,FC(R)))) = true;
    end
    Choose = Choose(size(P,1)+1:end);
end