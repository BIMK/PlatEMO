function [Xic,Xre] = DecisionVariableAnalysis(Problem,NCA,NIA)
% Decision variable analysis
% This code is modified from ControlVariableAnalysis.m and
% DividingDistanceVariables.m in MOEA/DVA

%------------------------------- Copyright --------------------------------
% Copyright 2017-2018 Yiping Liu
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------
   
    Xic = false(1,Problem.D);
    Xre = false(1,Problem.D);
    
    %% Find convergence-related variable
    for i = 1 : Problem.D
        x      = rand(1,Problem.D).*(Problem.upper-Problem.lower) + Problem.lower;
        S      = repmat(x,NCA,1);
        S(:,i) = ((1:NCA)'-1+rand(NCA,1))/NCA*(Problem.upper(i)-Problem.lower(i)) + Problem.lower(i);
        S = Problem.CalObj(S);
        S = unique(S,'rows'); % Delete the duplicate
        [~,MaxFNo] = NDSort(S,inf);
        if MaxFNo == size(S,1)
            Xic(i) = true;
        else
            Xre(i) = true;
        end
    end
    
    %% Interdependence 
    % Generate the initial population        
    PopDec = rand(Problem.N,Problem.D);
    PopDec = PopDec.*repmat(Problem.upper-Problem.lower,Problem.N,1) + repmat(Problem.lower,Problem.N,1);   
    PopObj = Problem.CalObj(PopDec);
    % Interdependence analysis
    interaction = false(Problem.D);
    interaction(logical(eye(Problem.D))) = true;
    for i = 1 : Problem.D-1
        for j = i+1 : Problem.D
            for time2try = 1 : NIA
                % Detect whether the i-th and j-th decision variables are interacting
                x    = randi(Problem.N);
                a2   = rand*(Problem.upper(i)-Problem.lower(i)) + Problem.lower(i);
                b2   = rand*(Problem.upper(j)-Problem.lower(j)) + Problem.lower(j);
                Decs = repmat(PopDec(x,:),3,1);
                Decs(1,i) = a2;
                Decs(2,j) = b2;
                Decs(3,[i,j]) = [a2,b2];
                F = Problem.CalObj(Decs);
                delta1 = F(1,:) - PopObj(x,:);
                delta2 = F(3,:) - F(2,:);
                interaction(i,j) = interaction(i,j) | any(delta1.*delta2<0);
                interaction(j,i) = interaction(i,j);                
            end
        end
    end
    
    %% Group based on Interdependence
    while sum(sum(interaction(Xic,Xre)))
        for i = find(Xic==1)
            fprintf('i=%d\n',i);
            if sum(interaction(i,Xre))
                Xic(i) = false;
                Xre(i) = true;
            end
        end      
    end
end