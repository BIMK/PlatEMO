function [Population,FrontNo,CrowdDis,result_sort] = Selection_ACRT(Population,N,a)
% Environmental selection based on the comprehensive ranking

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhiqiang Zeng (email: zhiqiang.zeng@outlook.com)

    [N1,~]       = size(Population.objs);  
    con_total    = sum(max(0,Population.cons),2);
    [FrontNo1,~] = NDSort([Population.objs,con_total],inf);
    CrowdDis1    = CrowdingDistance(Population.objs,FrontNo1);
    [~,r1]       = sortrows([FrontNo1',-CrowdDis1']);
    Rp(r1)       = 1 : N1;  % obtain ranking  based on CDP


    [FrontNo2,~] = NDSort(Population.objs,inf);
    CrowdDis2    = CrowdingDistance(Population.objs,FrontNo2);
    [~,r2]       = sortrows([FrontNo2',-CrowdDis2']);
    Rp1(r2)      = 1 : N1;  % obtain ranking  based on non-dominated sorting

    [~,randk1]   = sort(con_total); % Sorting based on the degree of constraint violation
    Rp2(randk1)  = 1 : N1;   
    R_sum        = a^2*Rp1+(1-a)^2*Rp2+(1-a)*Rp;

	[~,Rank] = sort(R_sum);
    result_sort(Rank) = 1 : N1;   
    
    Population = Population(Rank(1:N));
    FrontNo    = FrontNo1(Rank(1:N));
    CrowdDis   = CrowdDis1(Rank(1:N));  
end