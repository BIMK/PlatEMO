function [PopObjX,OffObjX,ArcObjX] = Distrubance(Problem,PopDec,OffDec,ArcDec)
% Re-evaluate solutions with the same disturbance

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N1 = size(PopDec,1);
    N2 = size(OffDec,1);
    P  = Problem.Perturb([PopDec;OffDec;ArcDec],1);
    PopObjX = P(1:N1).objs;
    OffObjX = P(N1+1:N1+N2).objs;
    ArcObjX = P(N1+N2+1:end).objs;
end