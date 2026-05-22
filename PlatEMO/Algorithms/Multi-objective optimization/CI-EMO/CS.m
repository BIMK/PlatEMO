function Popreal = CS(PopDec,PopObj,TSDec,TSObj)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Huixiang Zhen (email: zhenhuixiang@cug.edu.cn)

     % PopDec: candidate solutions, TSDec: modeling data
     % Non-dominated sorting of archived data
     [FrontNo,~] = NDSort(TSObj,inf);
     ND_TSObj    = TSObj(FrontNo==1,:);
     ND_TSDec    = TSDec(FrontNo==1,:);
     
     % Non-dominated Sorting of Candidate Solution Data
     [FrontNo1,~] = NDSort(PopObj,inf);
     ND_PopObj    = PopObj(FrontNo1==1,:);
     ND_PopDec    = PopDec(FrontNo1==1,:);   
     
     % parameters
     Popreal = [];
     rand1   = 1;

     %% Sampling Strategy
     if rand < rand1
         if ~isempty(PopDec)
             [Popreal_Dec1, Popreal_Obj1, index] = CI(ND_TSDec,ND_TSObj,ND_PopDec,ND_PopObj,TSDec,TSObj); % Composite indicator sampling
             Popreal = [Popreal; Popreal_Dec1];
         end
     end
end