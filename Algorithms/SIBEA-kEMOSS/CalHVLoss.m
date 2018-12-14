function HVLoss = CalHVLoss(PopObj,FrontNo)
% Calculate the weighted hypervolume (WHV) loss of each solution

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Liangli Zhen
  
    %% Calculate the WHV loss of each solution front by front
    HVLoss   = zeros(1,size(PopObj,1));
    RefPoint = max(PopObj,[],1) + 0.1;
    for f = setdiff(unique(FrontNo),inf)
        current  = find(FrontNo==f);
        totalWHV = CalHV(PopObj(current,:),RefPoint);
        for i = 1 : length(current)
            drawnow();
            currenti          = current([1:i-1,i+1:end]);
            HVLoss(current(i))= totalWHV - CalHV(PopObj(currenti,:),RefPoint);
        end
    end
end