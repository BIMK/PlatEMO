function [App,Dis] = subFitness(PopObj,P,Center,R)
% Update intersection point solution in each subregion

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    K          = length(R);
    [N,M]      = size(PopObj);

    % Normalize the population
    fmin       = min(PopObj,[],1);
    fmax       = max(PopObj,[],1);
    Obj        = (PopObj-repmat(fmin,N,1))./repmat(fmax-fmin,N,1);
    
    %% Calculate intersection point in each subregion
    if K == 1
        InterPoint = interPoint(Obj,P);
    else
        InterPoint = ones(N,M);
        
        % Allocation
        transformation = Allocation(Obj,Center,R);
       
        for i = 1 : K
            current     = find(transformation == i);
            if ~isempty(current)
                sInterPoint = interPoint(Obj(current,:),P(i,:));
                InterPoint(current,:) = sInterPoint;
            end
        end     
    end
    
    % Calculate the diversity and convergence of intersection points
    App = min(InterPoint-Obj,[],2);
%     App = sqrt(sum(InterPoint.^2,2)) - sqrt(sum(Obj.^2,2));

    % Calculate the diversity of each solution 
%     Dis = pdist2(InterPoint,InterPoint);
    Dis = distMax(InterPoint);
    Dis(logical(eye(length(Dis)))) = inf; 
end

function InterPoint = interPoint(PopObj,P)
% Calcualte the approximation degree of each solution, and the distances
% between the intersection points of the solutions

    [N,~] = size(PopObj);
    
    %% Calculate the intersections by gradient descent
    P     = repmat(P,N,1);      % Powers
    r     = ones(N,1);          % Parameters to be optimized
    lamda = zeros(N,1) + 0.002;   % Learning rates
    E     = sum((r.*PopObj).^P,2) - 1;   % errors
    for i = 1 : 1000
        newr = r - lamda.*E.*sum(P.*PopObj.^P.*r.^(P-1),2);
        newE = sum((newr.*PopObj).^P,2) - 1;
        update         = newr > 0 &sum(newE.^2) < sum(E.^2);
        r(update)      = newr(update);
        E(update)      = newE(update);
        lamda(update)  = lamda(update)*1.002; 
        lamda(~update) = lamda(~update)/1.002;
    end
    InterPoint = PopObj.*r;
end

function Dis = distMax(X)
% distMax pairwise distance between one set of observations
% Dis = distMax(X) returns a matrix D containing the maximum absolute
%   distance per dimension between each pair of observations in the MX-by-N
%   data matrix X and MX-by-N data matrix X. 

%   Example:
%      X = randn(100, 5);
%      D = distMax(X,Y);
%   >>size(D) = 100*100

    if isempty(X)
        error('X must be a non-empty matrix');
    end

    [N,~] = size(X); % nx,p
    Dis = zeros(N,N);
    for i = 1 : N
        for j = i+1 : N
            Dis(i,j) = max(abs(X(i,:)-X(j,:)));
        end
    end
    Dis = Dis + Dis';
end