function AngleM = AssociateAngle(PopObj)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [N,M]  = size(PopObj);
    cosine = 1 - pdist2(PopObj,PopObj,'cosine');
    cosine(logical(eye(length(cosine)))) = 0;
    Angle1 = acos(cosine);
    Angle2 = sort(Angle1,2);
    AngleM = zeros(N,1);
    Angle  = Angle2(:,1:M);
    for i = 1 : N
        for j = 1 : M
            AngleM(i) = AngleM(i) + Angle(i,j)/j;
        end
    end
end