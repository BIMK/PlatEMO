function d = Point2LineDistance(P, A, B)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Annibale Panichella

    d = zeros(size(P,1),1);
    for i = 1:size(P,1)
        pa = P(i,:) - A;
        ba = B - A;
        t = dot(pa, ba)/dot(ba, ba);
        d(i,1) = norm(pa - t * ba,2);
    end
end 