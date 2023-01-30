function [Population,FrontNo,CrowdDis] = EnvironmentalSelection(Population,N,a)
% Environmental selection

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
     
	[N1,~]       = size(Population.objs);   
    [FrontNo1,~] = NDSort(Population.objs,Population.cons,inf);
    
    CrowdDis1 = CrowdingDistance(Population.objs,FrontNo1);
  
    [~,r1] = sortrows([FrontNo1',-CrowdDis1']);
    Rc(r1) = 1 : N1;   %基于CDP获得的ranking

   [FrontNo2,~] = NDSort(Population.objs,0,inf);
    
    CrowdDis2 = CrowdingDistance(Population.objs,FrontNo2);
    
   [~,r2]  = sortrows([FrontNo2',-CrowdDis2']);
    Rp(r2) = 1 : N1;   %基于非支配排序获得的ranking
    
	R_sum = (1-a)*Rc+a*Rp;
    
	[~,Rank] = sort(R_sum);

    Population = Population(Rank(1:N));
    FrontNo    = FrontNo1(Rank(1:N));
    CrowdDis   = CrowdDis1(Rank(1:N));  
end