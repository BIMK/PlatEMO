function [App,Dis] = CalFitness(PopObj,P,A)
% Calcualte the approximation degree of each solution, and the distances
% between the intersection points of the solutions

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [N,M] = size(PopObj);

    %% Calculate the intersections by gradient descent
    A     = repmat(A,N,1);      % Coefficients
    P     = repmat(P,N,1);      % Powers
    r     = ones(N,1);          % Parameters to be optimized
    lamda = zeros(N,1) + 0.1;   % Learning rates
    E     = sum(A.*(PopObj.*repmat(r,1,M)).^P,2) - 1;   % errors
    for i = 1 : 1000
        newr = r - lamda.*E.*sum(A.*P.*PopObj.^P.*repmat(r,1,M).^(P-1),2);
        newE = sum(A.*(PopObj.*repmat(newr,1,M)).^P,2) - 1;
        update         = newr > 0 & sum(newE.^2) < sum(E.^2);
        r(update)      = newr(update);
        E(update)      = newE(update);
        lamda(update)  = lamda(update)*1.1;
        lamda(~update) = lamda(~update)/1.1;
    end
    PopObj1 = PopObj.*r;
    
    %% Calculate the convergence of each solution
    App = sqrt(sum(PopObj1.^2,2)) - sqrt(sum(PopObj.^2,2));

    %% Calculate the diversity of each solution
    Dis = pdist2(PopObj1,PopObj1);
    Dis(logical(eye(length(Dis)))) = inf;
end