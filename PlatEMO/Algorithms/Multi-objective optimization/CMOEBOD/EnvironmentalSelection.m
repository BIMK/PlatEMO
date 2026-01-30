function Next = EnvironmentalSelection(P,V)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhiyao Zhang (email: zhiyao.zhang.cn@gmail.com)

    %% Data Preprocess
    NV     = size(V,1);
    NP     = size(P.objs,1);
    PopObj = P.objs;                                                                                                                                                                                                                                            
    zmin   = min(PopObj,[],1); zmax = max(PopObj,[],1);
    PopObj = (PopObj - zmin)./max(zmax - zmin,10e-10);
    ObjMSE = P.objmse./(max(zmax - zmin,10e-10).^2);
    PopCon = P.cons; 
    ConMSE = P.conmse;

    %% Associate each solution to a reference vector
    Angle  = acos(1-pdist2(PopObj,V,'cosine'));
    Pindex = true(1,NP);
    Vindex = true(1,NV);
    
    %% Select one solution for each reference vector
    while any(Vindex)
        [~,associate] = min(Angle(Pindex,Vindex),[],2);
        Pexist = find(Pindex==1); 
        Vexist = find(Vindex==1); 
        for i = unique(associate)'
            current = find(associate==i);
            if isscalar(current)
                best = 1;
            else
                best = CPoB(PopObj(Pexist(current),:),ObjMSE(Pexist(current),:),...
                    PopCon(Pexist(current),:),ConMSE(Pexist(current),:),V(Vexist(i),:));
            end
            Pindex(Pexist(current(best))) = 0;
            Vindex(Vexist(i)) = 0;
        end
    end
    Next = Pindex==0;
end

function best = CPoB(PopObj,ObjMSE,PopCon,ConMSE,lamda)
    % Approximate Scalar
    PoF    = Feasible_Probability(PopCon,ConMSE);
    [N,M]  = size(PopObj);
    u      = lamda.*PopObj;
    sigma2 = (lamda.^2).*ObjMSE;
    sigma2 = abs(real(sigma2));
    mu     = u(:,1:2); sig2 = sigma2(:,1:2);
    [y,x]  = GPcal(mu,abs(sig2));
    if M >= 3
        for i = 3 : M
            mu = [y,u(:,i)]; sig2 = [x,sigma2(:,i)];
            [y,x] = GPcal(mu,abs(sig2));
        end
    end  
    y1 = sum(lamda.*PopObj,2);
    x1 = sum((lamda.^2).*ObjMSE,2);
    y  = y + 0.05*y1;
    x  = x + 0.05^2*x1;
   
    % Probabilistic Sorting
    sigma    = sqrt(x(reshape(ones(N,1)*(1:N),N*N,1),:) + repmat(x,N,1));
    mean     = y(reshape(ones(N,1)*(1:N),N*N,1),:) - repmat(y,N,1);
    mean(all(mean==0,2),:) = ones(sum(all(mean==0,2)),1);
    x_PD     = normcdf((0-mean)./abs(real(sigma)));
    y_PD     = 1 - x_PD;
    x_PD     = - x_PD.*PoF(reshape(ones(N,1)*(1:N),N*N,1),:);
    y_PD     = - y_PD.*repmat(PoF,N,1);
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
    dominate = sum(dominate,1);
    best     = find(dominate==min(dominate));
    if length(best) > 1
       best = best(randperm(length(best),1));
    end
end

function [y,x] = GPcal(mu,sig2)
    % Calculate the mu (x) and sigma^2 (y) of the aggregation function
    tao   = sqrt(sig2(:,1)+sig2(:,2));
    alpha = (mu(:,1)-mu(:,2))./tao;
    y     = mu(:,1).*normcdf(alpha) + mu(:,2).*normcdf(-alpha) + tao.*normpdf(alpha);
    x     = (mu(:,1).^2+sig2(:,1)).*normcdf(alpha) + (mu(:,2).^2+sig2(:,2)).*normcdf(-alpha)...
        + sum(mu,2).*tao.*normpdf(alpha) - y.^2;
end

function PoF = Feasible_Probability(PopCon,ConMSE)
    [N,M] = size(PopCon);
    PoF   = ones(N,1);
    for i = 1 : N
        for j = 1 : M
            PoF(i) = PoF(i) * normcdf((0-PopCon(i,j))/sqrt(ConMSE(i,j)));
        end
    end
end