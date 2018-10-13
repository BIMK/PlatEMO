function W = WeightingOptimization(x,t2,gamma,Lower,Upper)
% Weighting optimization of WOF based on ordered grouping and 0.2-value
% transformation

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Ordered grouping
    [~,rank] = sort(abs(x));
    [~,rank] = sort(rank);
    Group    = ceil(rank./ceil(length(x)./gamma));

    %% Generate random population
    N = 10;
    W     = rand(N,gamma)*2;
    x     = repmat(x,size(W,1),1);
    range = repmat(Upper-Lower,size(W,1),1);
    Population = INDIVIDUAL(x+0.2.*range.*(W(:,Group)-1));
    [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,N);
    
    %% Optimization
    for i = 2 : ceil(t2/N)
        MatingPool = TournamentSelection(2,N,FrontNo,-CrowdDis);
        OffspringW = Operator(W(MatingPool,:),0,2);
        Offspring  = INDIVIDUAL(x+0.2.*range.*(OffspringW(:,Group)-1));
        [Population,FrontNo,CrowdDis,Next] = EnvironmentalSelection([Population,Offspring],N);
        W = [W;OffspringW];
        W = W(Next,:);
    end
    [~,rank] = sortrows([FrontNo',-CrowdDis']);
    W        = W(rank(1),Group);
end