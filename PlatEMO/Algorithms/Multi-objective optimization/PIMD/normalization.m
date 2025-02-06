function [NormPop,NormComPop] = normalization(varargin)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yang Li (email: liyangwust@163.com)

    flag = 0;
    if nargin == 1
        PopObj = varargin{1};
        z = min(PopObj,[],1);
        if flag
            k = 1.5;
            [N,M] = size(PopObj);
            temp  = sort(PopObj);
            Q1    = temp(max(fix(N*0.25)),:);
            Q3    = temp(max(fix(N*0.75)),:);
            Max   = Q3+k*(Q3-Q1);
            znad = zeros(1,M);
            for i = 1 : M
                znad(i) = max(temp(temp(:,i)<Max(i),i));
            end
        else
            znad = max(PopObj,[],1);
        end
        % Normalization
        NormPop    = (PopObj-repmat(z,size(PopObj,1),1))./repmat(znad-z ,size(PopObj,1),1);
        NormComPop = [];
    else
        PopObj     = varargin{1};
        ComPopObj  = varargin{2};
        CombinePop = [PopObj;ComPopObj];

        z = min(CombinePop,[],1);
        if flag
            k = 1.5;
            [N,M] = size(CombinePop);
            temp  = sort(CombinePop);
            Q1    = temp(max(fix(N*0.25)),:);
            Q3    = temp(max(fix(N*0.75)),:);
            Max   = Q3+k*(Q3-Q1);
            znad = zeros(1,M);
            for i = 1 : M
                znad(i) = max(temp(temp(:,i)<Max(i),i));
            end
        else
            znad = max([PopObj;ComPopObj],[],1);
        end
        NormPop    = (PopObj-repmat(z,size(PopObj,1),1))./repmat(znad-z ,size(PopObj,1),1);
        NormComPop = (ComPopObj-repmat(z,size(ComPopObj,1),1))./repmat(znad-z ,size(ComPopObj,1),1);
    end
end