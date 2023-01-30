function P = subGFM(PopObj,Center,R,FrontNo)
% PF modeling for each subregion 

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    K     = size(Center,1);
    [N,M] = size(PopObj);

    % Normalize the population
    fmin = min(PopObj,[],1);
    fmax = max(PopObj,[],1);
    Obj  = (PopObj-repmat(fmin,N,1))./repmat(fmax-fmin,N,1);

    
    % PF modeling
    if K == 1
        P = GFM(Obj(FrontNo==1,:));
    else
        P = ones(K,M);
        
        % Allocation
        transformation = Allocation(Obj,Center,R);
        subFirstFront  = false(N,1);
        
        % Non-dominated sorting od each subregion
        for i = 1 : K
            current = find(transformation == i);
            if  ~isempty(current)
                [FNo,MFNo] = NDSort(PopObj(current,:),length(current));
                subFirstFront(current(FNo<MFNo|FNo==1)) = true;
            end
        end
        
        FTransformation = transformation(subFirstFront);
        PopObj          = PopObj(subFirstFront,:);
        RemainObj       = Obj(subFirstFront,:);        
        
        % GFM of each subregion
        if size(PopObj,1) > M
           for i = 1 : K
               current = find(FTransformation==i);
               if ~isempty(current)
                    if length(current) < M+1
                        [~,sDis] = sort(pdist2(RemainObj ,Center(i,:)));
                        current  = sDis(1:M+1);
                    end
                    p  = GFM(Obj(current,:));
                    P(i,:) = p;
               end
           end
        end
        
    end
    
end

function P = GFM(X)
% Generic front modeling

    [N,M] = size(X);
    X     = max(X,1e-12);
    P     = ones(1,M);
    lamda = 1;
    E     = sum(X.^repmat(P,N,1),2) - 1;
    MSE   = mean(E.^2);
    for epoch = 1 : 1000
        % Calculate the Jacobian matrix
        J = X.^repmat(P,N,1).*log(X);
        % Update the value of each weight
        while true
            Delta  = -(J'*J+lamda*eye(size(J,2)))^-1*J'*E;
            newP   = P + Delta(1:end)';
            newE   = sum(X.^repmat(newP,N,1),2) - 1;
            newMSE = mean(newE.^2);
            if newMSE < MSE && all(newP>1e-3)
                P     = newP;
                E     = newE;
                MSE    = newMSE;
                lamda = lamda/1.08;
                break;
            elseif lamda > 1e8
                return;
            else
                lamda = lamda*1.08;
            end
        end
    end
end