function phi = likehood(PDec,PCheby,theta,F) 
% The likehood calculation according to given theta values.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He

    N = size(PDec,1);
    R = zeros(N);
    for i = 1 : N
        for j = 1 : N
            R(i,j) = exp(-sum(theta.*(PDec(i,:)-PDec(j,:)).^2));
        end
    end
    bt     = (F'/R*F)\F'/R*PCheby;
    sigma2 = 1/N.*(PCheby-F*bt)'/R*(PCheby-F*bt);
    phi    = N.*log(sigma2)+log(det(R));
end