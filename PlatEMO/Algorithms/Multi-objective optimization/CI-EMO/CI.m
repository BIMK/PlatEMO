function [Popreal_Dec1, Popreal_Obj1, index] = CI(ND_TSDec,ND_TSObj,ND_PopDec,ND_PopObj,TSDec,TSObj)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Huixiang Zhen (email: zhenhuixiang@cug.edu.cn)

    % I1 
    angle     = acos(1-pdist2((max(ND_PopObj,0)),max(ND_TSObj,0),'cosine'));   
    [Angle,~] = min(angle,[],2);                                            
        
    % Objective Normalization
    ND_PopObj = (ND_PopObj - min(TSObj,[],1))./(max(TSObj,[],1)-min(TSObj,[],1));
    TSObj     = (TSObj- min(TSObj,[],1))./(max(TSObj,[],1)-min(TSObj,[],1));
    Z         = min(TSObj,[],1);

    % I2 
    ddt         = pdist2(ND_PopObj, TSObj,'euclidean');
    ddt(ddt==0) = inf;
    CN          = min(ddt,[],2);

    % I3
    DC = pdist2(ND_PopObj, Z,'euclidean');

    % Indicators Normalization
    N_angle = (Angle - min(Angle,[],1))./(max(Angle,[],1)-min(Angle,[],1));
    N_CN    = (CN - min(CN,[],1))./(max(CN,[],1)-min(CN,[],1));
    N_DC    = (DC - min(DC,[],1))./(max(DC,[],1)-min(DC,[],1));
    
    % CI
    r1 = rand;
    r2 = rand;
    r3 = rand;
    CI = r1 * N_CN + r2 * N_angle - r3 * N_DC;
    
    % Selection
    [~, index]   = max(CI);
    Popreal_Dec1 = ND_PopDec(index,:);
    Popreal_Obj1 = ND_PopObj(index,:);
end