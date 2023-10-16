function optimo = LocalSearch(Problem,pos)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    MaxIter = 5;
    Tol     = 1e-3;
    step    = 1;
    k       = 1;
    error   = 10;
    while error>Tol && k<MaxIter
        grad(1,:)    = FiniteDifference(pos,Problem);
        offspringdec = pos.dec - step*grad(1,:);
        offspringdec = min(max(offspringdec,Problem.lower),Problem.upper);
        offspring    = Problem.Evaluation(offspringdec);
        grad(2,:)    = FiniteDifference(offspring,Problem);
        step         = abs((offspring.dec-pos.dec)*(grad(2,:)-grad(1,:))')/norm(grad(2,:)-grad(1,:))^2;
        error = norm(offspring.dec-pos.dec);
        pos   = offspring;
        k     = k + 1;
    end
    optimo = pos;
end

function df = FiniteDifference(X,Problem)
    if any(X.con>0)
        df = Problem.CalConGrad(X.dec)';
        df = sum(df,2);
    else
        df = Problem.CalObjGrad(X.dec)';
    end
end