function Distance = CalDistance(PopObj,RefPoint)
% Calculate the distance between each solution to each adjusted reference
% point

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N  = size(PopObj,1);
    NR = size(RefPoint,1);
    PopObj   = max(PopObj,1e-6);
    RefPoint = max(RefPoint,1e-6);

    %% Adjust the location of each reference point
    Cosine = 1 - pdist2(PopObj,RefPoint,'cosine');
    NormR  = sqrt(sum(RefPoint.^2,2));
    NormP  = sqrt(sum(PopObj.^2,2));
    d1     = repmat(NormP,1,NR).*Cosine;
    d2     = repmat(NormP,1,NR).*sqrt(1-Cosine.^2);
    [~,nearest] = min(d2,[],1);
    RefPoint    = RefPoint.*repmat(d1(N.*(0:NR-1)+nearest)'./NormR,1,size(RefPoint,2));
    
    %% Calculate the distance between each solution to each point
    Distance = pdist2(PopObj,RefPoint);
end