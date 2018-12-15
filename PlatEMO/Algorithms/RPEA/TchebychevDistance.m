function Distance = TchebychevDistance(PopObj,R)
% Calculate the Tchebychev distance between each individual and reference
% point

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    fmax     = max(PopObj,[],1);
    fmin     = min(PopObj,[],1);
    Distance = zeros(size(PopObj,1),size(R,1));
    for i = 1 : size(PopObj,1)
        Distance(i,:) = max((repmat(PopObj(i,:),size(R,1),1)-R)./repmat(fmax-fmin,size(R,1),1)./size(PopObj,2),[],2)';
    end
end