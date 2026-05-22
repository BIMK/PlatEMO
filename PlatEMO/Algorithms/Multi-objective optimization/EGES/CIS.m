function Popreal = CIS(PopDec,PopObj,TSDec,TSObj)
% PopDec: candidate solution data; TSDec: archived data

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Huixiang Zhen (email: zhenhuixiang@cug.edu.cn)

     %% Non-dominated sorting of archived data
     [FrontNo,~] = NDSort(TSObj,inf);
     ND_TSObj    = TSObj(FrontNo==1,:);
     ND_TSDec    = TSDec(FrontNo==1,:);
     
     %% Non-dominated sorting of candidate solution data
     [FrontNo1,~] = NDSort(PopObj,inf);
     ND_PopObj    = PopObj(FrontNo1==1,:);
     ND_PopDec    = PopDec(FrontNo1==1,:);   
     
     %% parameters
     Popreal = [];
     rand1   = 1;

     %% CI sampling
     if rand < rand1
         if ~isempty(PopDec)
             CI_predicted = CI(ND_PopObj,TSObj);
             [~, index]   = max(CI_predicted);
             Popreal      = [Popreal; PopDec(index,:)];
         end
     end
end


