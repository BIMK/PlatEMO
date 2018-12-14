function FC = CalFC(PopObj,Zmin)
% Calculate the FC value of each solution

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [N,M]  = size(PopObj);
    PopObj = PopObj - repmat(Zmin,N,1);
    
    %% Favorable weight
    w     = zeros(N,M);
    bound = any(PopObj==repmat(Zmin,N,1),2);
    w(repmat(bound,1,M) & PopObj==0) = 1;
    w(repmat(bound,1,M) & PopObj~=0) = 0;
    w(~bound,:) = 1./PopObj(~bound,:)./repmat(sum(1./PopObj(~bound,:),2),1,M);
    
    %% Calculate the FC value
    FC = max(w.*PopObj,[],2);
    FC = max(FC,1e-6);
end