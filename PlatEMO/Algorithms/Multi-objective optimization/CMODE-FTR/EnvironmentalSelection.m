function [Population] = EnvironmentalSelection(Population,N,a,Problem)
% Environmental selection

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhiqiang Zeng (email: zhiqiang.zeng@outlook.com)

    [N1,~]       = size(Population.objs);
    [FrontNo1,~] = NDSort(Population.objs,Population.cons,inf);
    CrowdDis1    = CrowdingDistance(Population.objs,FrontNo1);
    [~,r1]       = sortrows([FrontNo1',-CrowdDis1']);
    Rc(r1)       = 1 : N1;   % Obtain the ranking of each individual based on CDP
    
    FrontNo2  = NDSort(Population.objs,0,inf);
    CrowdDis2 = CrowdingDistance(Population.objs,FrontNo2);
    [~,r2]    = sortrows([FrontNo2',-CrowdDis2']);
    Rp(r2)    = 1 : N1;   % Obtain the ranking of each individual based on non-dominated sorting
    
    
    pro_l = 1-length(find(sum(max(0,Population.cons),2)>0))/N1; % Calculate feasibility rate
    b     = 1.0/((Problem.FE/Problem.maxFE)^2+1)-0.5;
    r     = a*(1-b)+pro_l*b;
    R_sum = (1-r)*Rc+r*Rp; % To fuse two rankings.
    
    [~,Rank]   = sort(R_sum);
    Population = Population(Rank(1:N));
end