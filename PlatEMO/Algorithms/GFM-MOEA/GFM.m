function [P,A] = GFM(X)
% Generic front modeling

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [N,M] = size(X);
    X     = max(X,1e-12);
    P     = ones(1,M);
    A     = ones(1,M);
    lamda = 1;
	E     = sum(repmat(A,N,1).*X.^repmat(P,N,1),2) - 1;
    MSE   = mean(E.^2);
    for epoch = 1 : 1000
        % Calculate the Jacobian matrix
        J = [repmat(A,N,1).*X.^repmat(P,N,1).*log(X),X.^repmat(P,N,1)];
        % Update the value of each weight
        while true
            Delta  = -(J'*J+lamda*eye(size(J,2)))^-1*J'*E;
            newP   = P + Delta(1:M)';
            newA   = A + Delta(M+1:end)';
            newE   = sum(repmat(newA,N,1).*X.^repmat(newP,N,1),2) - 1;
            newMSE = mean(newE.^2);
            if newMSE < MSE && all(newP>1e-3) && all(newA>1e-3)
                P     = newP;
                A     = newA;
                E     = newE;
                MSE   = newMSE;
                lamda = lamda/1.1;
                break;
            elseif lamda > 1e8
                return;
            else
                lamda = lamda*1.1;
            end
        end
    end
end