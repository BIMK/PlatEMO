function [WMM,VADS] = CalMetric(PopObj,W)
% Calculate the weighted min-max metric and the vector angle distance
% scaling metric between each solution and vector

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

	N  = size(PopObj,1);
    NW = size(W,1);

    %% Normalization
    Zmax   = max(PopObj,[],1);
    Zmin   = min(PopObj,[],1);
    PopObj = (PopObj-repmat(Zmin,N,1))./(repmat(Zmax-Zmin,N,1));
    W      = 1./(W-repmat(Zmin,NW,1)+eps);
    
    %% Calculate WMM
    WMM = zeros(N,NW);
    for i = 1 : NW
        WMM(:,i)  = max(PopObj.*repmat(W(i,:)./norm(W(i,:)),N,1),[],2);
    end
    
    %% Calculate VADS
    if nargout > 1
        VADS  = zeros(N,NW);
        NormP = sqrt(sum(PopObj.^2,2));
        for i = 1 : NW
            VADS(:,i) = NormP./(PopObj*W(i,:)'./NormP).^100;
        end
    end
end