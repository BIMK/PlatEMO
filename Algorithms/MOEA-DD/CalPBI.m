function PBI = CalPBI(PopObj,W,Region,Z,Sub)
% Calculate the PBI value between each solution and its associated weight
% vector

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Z   = repmat(Z,sum(Sub),1);
    NormW = sqrt(sum(W(Region(Sub),:).^2,2));
    d1  = abs(sum((PopObj(Sub,:)-Z).*W(Region(Sub),:),2))./NormW;
    d2  = sqrt(sum((PopObj(Sub,:)-(Z+W(Region(Sub),:).*repmat(d1./NormW,1,size(W,2)))).^2,2));
    PBI = zeros(1,size(PopObj,1));
    PBI(Sub) = d1 + 5*d2;
end