function [Population,Groups,K] = spiltVariables(Problem,BU,BD,s)
% Decision variable grouping with RDG2 in SACCEAMII

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Initialize
    D = length(BU);
    Groups     = zeros(1,D);
    Population = [];
    Separable  = [];
    NonSeparable = [];
    
    %% Recursive Decomposition Grouping Methods
    xLL  = BD;
    xNew = Problem.Evaluation(xLL);
    yLL  = xNew.objs;
    Population = [Population,xNew];
    X1   = 1;
    X2   = 2:D;
    % Process
    while ~isempty(X2)
        % Interaction detection
        [sub1,Population] = INTERACT(Problem,Population,X1,X2,xLL,yLL,BU,BD);
        if numel(sub1) == numel(X1) % isequal(sub1,X1)
            if numel(X1) == 1
                Separable = [Separable,X1];
            else
                NonSeparable = [NonSeparable,X1];
            end
			X1 = X2(1);
			X2 = X2(2:end);
        else
            X1 = sub1;
            X2 = X2(~ismember(X2,X1));
        end
        
        if isempty(X2)
            if numel(X1) > 1
                NonSeparable = [NonSeparable,X1];
            else
                Separable = [Separable,X1];
            end
        end
    end
    
    %% Split the separable decision varables
     if ~isempty(Separable)
         num = floor(numel(Separable)/s);
         if num <= 1
             K   = 1;
             Groups(Separable) = K;
         else
             K = num;
             idx = 1:1:numel(Separable);
             for i = 1 : num-1
                 select = randperm(numel(idx),s);
                 Groups(Separable(idx(select))) = i;
                 idx(select) = [];
             end
             Groups(Separable(idx)) = num;
         end
         if ~isempty(NonSeparable)
             K = K + 1;
             Groups(NonSeparable) = K;
         end
     else
         K = 1;
         Groups(NonSeparable) = K;
     end
end

function [X1,Population] = INTERACT(Problem,Population,X1,X2,xLL,yLL,BU,BD)
    % Determine gamma
    muM   = eps/2;
    gamma = @(n)(n.*muM)./(1-n.*muM);
    
    xUL     = xLL;
    xUL(X1) = BU(X1);
    
    % Calculate delta1
    tNew       = Problem.Evaluation(xUL);
    delta1     = yLL - tNew.objs;
    Population = [Population,tNew];
    
    % Calculate delta2
    xLM     = xLL;
    xLM(X2) = (BU(X2)+BD(X2))/2;
    xUM     = xUL;
    xUM(X2) = (BU(X2)+BD(X2))/2;
    tNew1   = Problem.Evaluation(xLM);
    tNew2   = Problem.Evaluation(xUM);
    delta2  = tNew1.objs - tNew2.objs;
    Population = [Population,tNew1,tNew2];
    % Update
    F4 = [yLL,tNew.objs,tNew1.objs,tNew2.objs];
    threshold = gamma(numel(BU).^0.5+2)*sum(abs(F4));
    if abs(delta1-delta2) > threshold
        if numel(X2) == 1
            X1 = union(X1,X2);
        else
            % Divide X2 into equally-sized groups G1, G2
            mid = floor(length(X2)/2);
            G1  = X2(1:mid);
            G2  = X2(mid+1:end);
            [subX1,Population] = INTERACT(Problem,Population,X1,G1,xLL,yLL,BU,BD);
            [subX2,Population] = INTERACT(Problem,Population,X1,G2,xLL,yLL,BU,BD);
            X1  = union(subX1,subX2);
        end
    end
end