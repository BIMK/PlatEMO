function [Subcomponents,Population] = DividingDistanceVariables(Problem,NIA,DiverIndexes,ConverIndexes)
% Dividing distance variables

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    %% Generate the initial population
    PopDec = zeros(Problem.N,Problem.D);
    % Generate the diverse variables by uniformly sampling method
    if sum(DiverIndexes) == 1
        PopDec(:,DiverIndexes) = (0:Problem.N-1)/(Problem.N-1);
    elseif sum(DiverIndexes) > 4
        PopDec(:,DiverIndexes) = rand(Problem.N,sum(DiverIndexes));
    else
        PopDec(:,DiverIndexes) = UDall(Problem.N,sum(DiverIndexes));
    end
    % Randomly generate the distance variables
    PopDec(:,ConverIndexes) = rand(Problem.N,sum(ConverIndexes));
    % Generate the initial population
    PopDec     = PopDec.*repmat(Problem.upper-Problem.lower,Problem.N,1) + repmat(Problem.lower,Problem.N,1);
    Population = Problem.Evaluation(PopDec);
    
    %% Interdependence analysis
    interaction = false(Problem.D);
    interaction(logical(eye(Problem.D))) = true;
    for i = 1 : Problem.D-1
        for j = i+1 : Problem.D
            drawnow('limitrate');
            for time2try = 1 : NIA
                % Detect whether the i-th and j-th decision variables are
                % interacting
                x    = randi(Problem.N);
                a2   = rand*(Problem.upper(i)-Problem.lower(i)) + Problem.lower(i);
                b2   = rand*(Problem.upper(j)-Problem.lower(j)) + Problem.lower(j);
                Decs = repmat(Population(x).dec,3,1);
                Decs(1,i) = a2;
                Decs(2,j) = b2;
                Decs(3,[i,j]) = [a2,b2];
                F = Problem.Evaluation(Decs);
                delta1 = F(1).obj - Population(x).obj;
                delta2 = F(3).obj - F(2).obj;
                interaction(i,j) = interaction(i,j) | any(delta1.*delta2<0);
                interaction(j,i) = interaction(i,j);
                % Update the solution
                if ConverIndexes(j) && all(F(2).obj <=Population(x).obj)
                    Population(x) = F(2);
                end
                if ConverIndexes(i) && all(F(1).obj <=Population(x).obj)
                    Population(x) = F(1);
                end
                if all(ConverIndexes([i,j])) && all(F(3).obj <=Population(x).obj)
                    Population(x) = F(3);
                end
            end
        end
    end

    %% Dividing distance variables
    Subcomponents = {};
    divided = false(1,Problem.D);
    while ~all(divided(ConverIndexes))
        x = find(~divided & ConverIndexes,1);
        while sum(any(interaction(x,ConverIndexes),1)) > length(x)
            x = find(any(interaction(x,:),1) & ConverIndexes);
        end
        Subcomponents = [Subcomponents,x];
        divided(x) = true;
    end
end