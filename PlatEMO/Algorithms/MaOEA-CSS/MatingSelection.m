function MatingPool = MatingSelection(PopObj,Zmin)
% The mating selection of MaOEA-CSS

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [N,M] = size(PopObj);
    
    %% Calculate the ASF value of each solution
    W      = max(1e-6,PopObj./repmat(sum(PopObj,2),1,M));
    PopObj = PopObj - repmat(Zmin,N,1);
    ASF    = max(PopObj./W,[],2);
    % Obtain the rank value of each solution's ASF value
    [~,rank]    = sort(ASF);
    [~,ASFRank] = sort(rank);
    
    %% Calculate the minimum angle of each solution to others
    Angle = acos(1-pdist2(PopObj,PopObj,'cosine'));
    Angle(logical(eye(N))) = inf;
    Amin  = min(Angle,[],2);
    
    %% Binary tournament selection
    MatingPool = zeros(1,N);
    for i = 1 : N
        p = randperm(N,2);
        if ASF(p(1)) < ASF(p(2)) && Amin(p(1)) > Amin(p(2))
            p = p(1);
        else
            p = p(2);
        end
        if rand < 1.0002+ASFRank(p)/N
            MatingPool(i) = p;
        else
            MatingPool(i) = randi(N);
        end
    end
end