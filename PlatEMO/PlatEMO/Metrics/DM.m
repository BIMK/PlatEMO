function Score = DM(PopObj,PF)
% <metric> <max>
% Metric for diversity

%------------------------------- Reference --------------------------------
% K. Deb and S. Jain, Running performance metrics for evolutionary
% multi-objective optimization, KanGAL Report 2002004, 2002.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    fmax  = max(PF,[],1);
    fmin  = min(PF,[],1);
    H     = calGrid(PF(:,1:end-1),fmax(1:end-1),fmin(1:end-1),size(PopObj,1));
    h     = H & calGrid(PopObj(:,1:end-1),fmax(1:end-1),fmin(1:end-1),size(PopObj,1));
    Score = calM(h,H)./calM(H,H);
end

function h = calGrid(P,fmax,fmin,div)
% Determine whether each grid has at least one point

    [N,M] = size(P);
    d     = (fmax-fmin)./div;
    GLoc  = ceil((P-repmat(fmin,N,1))./repmat(d,N,1));
    GLoc  = max(1,GLoc);
    h     = zeros(M,div);
    for i = 1 : M
        h(i,:) = ismember(1:div,GLoc(:,i));
    end
end

function m = calM(h,H)
% Calculate the value function m()

    M = size(h,1);
    h = [ones(M,1),h,ones(M,1)];
    H = [ones(M,1),H,ones(M,1)];
    m = 0;
    for i = 1 : M
        for j = 2 : size(h,2)-1
            if H(i,j)
                if h(i,j)
                    if h(i,j-1)
                        if h(i,j+1)
                            m = m + 1;
                        else
                            m = m + 0.67;
                        end
                    else
                        if h(i,j+1)
                            m = m + 0.67;
                        else
                            m = m + 0.75;
                        end
                    end
                else
                    if h(i,j-1)
                        if h(i,j+1)
                            m = m + 0.75;
                        else
                            m = m + 0.5;
                        end
                    else
                        if h(i,j+1)
                            m = m + 0.5;
                        else
                            m = m + 0;
                        end
                    end
                end
            end
        end
    end
end