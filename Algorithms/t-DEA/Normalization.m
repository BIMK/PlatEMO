function [PopObj,z,znad] = Normalization(PopObj,z,znad)
% Normalize the population and update the ideal point and the nadir point

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [N,M] = size(PopObj);

    %% Update the ideal point
    z = min(z,min(PopObj,[],1));
    
    %% Update the nadir point
    % Identify the extreme points
    W = zeros(M) + 1e-6;
    W(logical(eye(M))) = 1;
    ASF = zeros(N,M);
    for i = 1 : M
        ASF(:,i) = max(abs((PopObj-repmat(z,N,1))./(repmat(znad-z,N,1)))./repmat(W(i,:),N,1),[],2);
    end
    [~,extreme] = min(ASF,[],1);
    % Calculate the intercepts
    Hyperplane = (PopObj(extreme,:)-repmat(z,M,1))\ones(M,1);
    a = (1./Hyperplane)' + z;
    if any(isnan(a)) || any(a<=z)
        a = max(PopObj,[],1);
    end
    znad = a;
    
    %% Normalize the population
    PopObj = (PopObj-repmat(z,N,1))./(repmat(znad-z,N,1));
end