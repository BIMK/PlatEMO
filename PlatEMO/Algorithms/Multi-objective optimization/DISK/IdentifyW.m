function [W,ideal] = IdentifyW(DB,N,M)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhiyao Zhang (email: z.zhang0@csu.edu.cn)

    %% Preparing Data 
    V      = UniformPoint(10*N,M);
    A2Obj  = DB.objs;
    [F_,~] = NDSort(A2Obj,1);
    A2Obj  = A2Obj(F_==1,:);

    %% Translate Coordinate
    nadir = max(A2Obj,[],1);
    ideal = min(A2Obj,[],1);
    ideal = ideal - (nadir-ideal)/10 - 0.1*ones(1,M);
    A2Obj = A2Obj - ideal; 
    
    %% Calculate Angle between A2Obj and V
    Angle  = acos(1-pdist2(V,A2Obj,'cosine'));
    Angle_ = sort(Angle,2);
    index  = find(Angle_(:,1)==max(Angle_(:,1)));
    
    %% Identify the Farthest W
    if length(index) > 1
        index = index(randperm(length(index),1));
    end
    W = V(index,:);  
end