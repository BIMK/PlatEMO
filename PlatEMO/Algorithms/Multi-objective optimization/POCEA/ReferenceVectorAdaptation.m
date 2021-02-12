function [V1,V2] = ReferenceVectorAdaptation(PopObj,V1,V2)
% Reference vector adaption strategy of two reference vector sets

% Copyright (c) 2020-2021 Cheng He

    PopObj = max(PopObj,[],1)-min(PopObj,[],1);
    V1 	   = V1.*repmat(PopObj,size(V1,1),1); 
    V2     = V2.*repmat(PopObj,size(V2,1),1); 
end