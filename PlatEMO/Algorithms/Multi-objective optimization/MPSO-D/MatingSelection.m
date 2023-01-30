function [Parent,Pbest,Gbest] = MatingSelection(PopObj,B,W,Z)
% Mating selection of MPSO/D

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N = size(PopObj,1);

    %% Select N parents according to their crowding distances
    PopObj = PopObj - repmat(Z,size(PopObj,1),1);
    Parent = TournamentSelection(2,N,-CrowdingDistance(PopObj));
    
    %% Select the pbest and gbest of each parent
    Pbest = zeros(size(Parent));
    Gbest = zeros(size(Parent));
    for i = 1 : N
        if rand < 0.9
            P = B(i,:);
        else
            P = 1 : N;
        end
        Pbest(i) = P(randi(length(P)));
        [~,best] = max(PopObj(P,:)*mean(W(P,:),1)'./sum(PopObj(P,:).^2,2).^0.6);
        Gbest(i) = P(best);
    end
end