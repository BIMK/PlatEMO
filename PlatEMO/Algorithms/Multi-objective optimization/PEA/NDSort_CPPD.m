function [FrontNo,MaxFNo] = NDSort_CPPD(PopObj,ObjMSE,PopCon,ConMSE,nSort)
% Constrained Probabilistic Dominance (CPD)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhiyao Zhang (email: zhiyao_zhang0@163.com)

    N     = size(PopObj,1);
    LFP   = Feasible_Probability(PopCon,ConMSE);
    sigma = sqrt(ObjMSE(reshape(ones(N,1)*(1:N),N*N,1),:) + repmat(ObjMSE,N,1));
    mean  = PopObj(reshape(ones(N,1)*(1:N),N*N,1),:) - repmat(PopObj,N,1);
    x_PD  = normcdf((0-mean)./sigma);
    y_PD  = 1 - x_PD;
    x_PD  = - x_PD.*LFP(reshape(ones(N,1)*(1:N),N*N,1),:);
    y_PD  = - y_PD.*repmat(LFP,N,1);
    
    dominate = false(N);
    for i = 1 : N-1
        for j = i+1 : N
            if all(x_PD(N*(i-1)+j,:) <= y_PD(N*(i-1)+j,:)) && ~all(x_PD(N*(i-1)+j,:) == y_PD(N*(i-1)+j,:))
                dominate(i,j) = true;
            elseif all(x_PD(N*(i-1)+j,:) >= y_PD(N*(i-1)+j,:)) && ~all(x_PD(N*(i-1)+j,:) == y_PD(N*(i-1)+j,:))
                dominate(j,i) = true;
            end
        end
    end

    FrontNo = inf(1,N);
    MaxFNo  = 0;
    while sum(FrontNo~=inf) < min(nSort,N)
        MaxFNo                     = MaxFNo + 1;
        current                    = find(FrontNo==inf);
        dominate_                  = sum(dominate(current,current),1);
        index                      = find(dominate_==min(dominate_));
        FrontNo(current(index))    = MaxFNo;
        dominate(current(index),:) = false;
    end
end

function LFP = Feasible_Probability(PopCon,ConMSE)
    [N,M] = size(PopCon);
    LFP   = ones(N,1);
    for i = 1 : N
        for j = 1 : M
             LFP(i) = min([LFP(i),normcdf((0-(PopCon(i,j)+0*sqrt(ConMSE(i,j))))/sqrt(ConMSE(i,j)))]);
        end
    end
end