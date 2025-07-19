function[FrontNo,MaxFNo] = NDSort_PDPD(varargin)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhiyao Zhang (email: zhiyao_zhang0@163.com)

    %% Read data
    PopObj = varargin{1};
    ObjMSE = varargin{2};
    if nargin == 3
        nSort = varargin{3};
        mark  = 1;
    elseif nargin == 5
        PopCon      = varargin{3};
        ConMSE      = varargin{4};
        nSort       = varargin{5};
        [LPoF,TPoF] = Feasible_Probability(PopCon,ConMSE);
        mark        = 0;
    end
    epsilon = 0.75;

    %% Obtain Dominance Matrix 
    N        = size(PopObj,1);
    sigma    = sqrt(ObjMSE(reshape(ones(N,1)*(1:N),N*N,1),:) + repmat(ObjMSE,N,1));
    mean     = PopObj(reshape(ones(N,1)*(1:N),N*N,1),:) - repmat(PopObj,N,1);
    x_PD     = normcdf((0-mean)./sigma);
    y_PD     = 1 - x_PD;
    dominate = false(N);
    for i = 1 : N-1
        for j = i+1 : N
            Pi             = x_PD(N*(i-1)+j,:);
            Pj             = y_PD(N*(i-1)+j,:);
            index1         = find(abs(Pi - Pj)<=epsilon);
            index2         = 1 : length(Pi);
            index2(index1) = [];
            PDi            = prod(Pi(index1));
            PDj            = prod(Pj(index1));
    
            if mark == 1
                % PDPD
                if all([-Pi(index2),-PDi] <= [-Pj(index2),-PDj]) && ~all([-Pi(index2),-PDi] == [-Pj(index2),-PDj])
                    flag = 1;
                elseif all([-Pi(index2),-PDi] >= [-Pj(index2),-PDj]) && ~all([-Pi(index2),-PDi] == [-Pj(index2),-PDj])
                    flag = 2;
                else
                    flag = 3;
                end
            else 
                % Constrained PDPD
                if ((LPoF(i) >= 0.5) && (LPoF(j) >= 0.5)) || (LPoF(i) == LPoF(j))
                    if all([-Pi(index2),-PDi] <= [-Pj(index2),-PDj]) && ~all([-Pi(index2),-PDi] == [-Pj(index2),-PDj])
                        flag = 1;
                    elseif all([-Pi(index2),-PDi] >= [-Pj(index2),-PDj]) && ~all([-Pi(index2),-PDi] == [-Pj(index2),-PDj])
                        flag = 2;
                    else
                        flag = 3;
                    end
                else
                    [~,flag] = max([TPoF(i),TPoF(j)]);
                end
            end
            
            if flag == 1
                dominate(i,j) = true;
            elseif flag == 2
                dominate(j,i) = true;
            end
        end
    end
    
    %% Sorting
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

function [LPoF,TPoF] = Feasible_Probability(PopCon,ConMSE)
    [N,M] = size(PopCon);
    LPoF  = ones(N,1);
    TPoF  = ones(N,1);
    for i = 1 : N
        for j = 1 : M
             LPoF(i) = min([LPoF(i),normcdf((0-PopCon(i,j))/sqrt(ConMSE(i,j)))]);
             TPoF(i) = TPoF(i)*normcdf((0-PopCon(i,j))/sqrt(ConMSE(i,j)));
        end
    end
end