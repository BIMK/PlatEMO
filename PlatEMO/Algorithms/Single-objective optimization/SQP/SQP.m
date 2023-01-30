classdef SQP < ALGORITHM
% <single> <real> <large/none> <constrained/none>
% Sequential quadratic programming

%------------------------------- Reference --------------------------------
% P. T. Boggs and J. W. Tolle, Sequential quadratic programming, Acta
% Numerica, 1995, 4(1): 1-51.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Generate random solution
            X     = Problem.Initialization(1);
            ro    = 0.5;
            eta   = 0.1;
            lam   = zeros(1,length(X.con));
            Bk    = eye(Problem.D);
            sigma = 0.8;
            dfk   = Problem.CalObjGrad(X.dec);
            Ai    = Problem.CalConGrad(X.dec);
            
            %% Optimization
            while Algorithm.NotTerminated(X)
                [~,dk,mu,~] = qpsubp(dfk',Bk,[],[],-Ai,-X.con');
                tau = max(norm(mu,inf),norm(lam,inf));
                if sigma*(tau+0.05) >= 1
                    sigma = 1/(tau+2*0.05);
                end
                for mk = 0 : 20
                    temp = eta*ro^mk*dphi1(X,sigma,dfk,dk);
                    X1   = Problem.Evaluation(X.dec+ro^mk*dk');
                    if phi1(X1,sigma)-phi1(X,sigma) < temp
                        break;
                    end
                end
                [dfk0,Ai0] = deal(dfk,Ai);
                dfk        = Problem.CalObjGrad(X1.dec);
                Ai         = Problem.CalConGrad(X1.dec);
                lam = pinv(-Ai')*dfk';
                sk  = (X1.dec-X.dec)';
                yk  = dlax(dfk,-Ai,lam) - dlax(dfk0,-Ai0,lam);
                if sk'*yk>0.2*sk'*Bk*sk
                    omega = 1;
                else
                    omega = 0.8*sk'*Bk*sk/(sk'*Bk*sk-sk'*yk);
                end
                zk = omega*yk + (1-omega)*Bk*sk;
                Bk = Bk + zk*zk'/(sk'*zk) - (Bk*sk)*(Bk*sk)'/(sk'*Bk*sk);
                X  = X1;
            end
        end
    end
end

function p = phi1(X,sigma)
    p = X.obj + 1/sigma*norm(max(X.con',0),1);
end

function dp = dphi1(X,sigma,df,d)
    dp = df*d - 1/sigma*norm(max(X.con',0),1);
end

function dl = dlax(df,Ai,lam)
    dl = df' - Ai'*lam;
end