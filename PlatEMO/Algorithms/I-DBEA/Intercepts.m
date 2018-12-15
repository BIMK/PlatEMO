function a = Intercepts(PopObj)
% Calculate the intercepts used in normalization

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [N,M] = size(PopObj);

    %% Find the extreme points
    [~,Choosed(1:M)] = min(PopObj,[],1);
    L2NormABO        = zeros(N,M);
    for i = 1 : M
    	L2NormABO(:,i) = sum(PopObj(:,[1:i-1,i+1:M]).^2,2);
    end
    [~,Choosed(M+1:2*M)] = min(L2NormABO,[],1);
    [~,Extreme]          = max(PopObj(Choosed,:),[],1);
    Extreme              = unique(Choosed(Extreme));
    
    %% Calculate the intercepts
    if length(Extreme) < M
        a = max(PopObj,[],1);
    else
        Hyperplane = PopObj(Extreme,:)\ones(M,1);
        a = 1./Hyperplane';
    end
end