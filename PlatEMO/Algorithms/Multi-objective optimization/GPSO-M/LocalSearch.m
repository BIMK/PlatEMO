function optimo = LocalSearch(Problem,pos,w)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
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
        grad(1,:)    = FiniteDifference(pos,w,Problem);
        offspringdec = pos.dec - step*grad(1,:);
        offspringdec = min(max(offspringdec,Problem.lower),Problem.upper);
        offspring    = Problem.Evaluation(offspringdec);
        grad(2,:)    = FiniteDifference(offspring,w,Problem);
        step         = abs((offspring.dec-pos.dec)*(grad(2,:)-grad(1,:))')/norm(grad(2,:)-grad(1,:))^2;
        error = norm(offspring.dec-pos.dec);
        pos   = offspring;
        k     = k + 1;
    end
    optimo = pos;
end

function df = FiniteDifference(X,W,Problem)
    feasible = X.con <= 0;
    if ~all(feasible)
        [~,df] = Problem.CalGrad(X.dec);
        df(feasible,:) = 0;
        df = sum(df',2);
    else
        df = Problem.CalGrad(X.dec)';
        df = df*W';
    end
end