function Fitness = CalFitness(PopObj,PopCon,priority,current_cons,constraint_handing)
% Calculate the fitness of each solution

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N = size(PopObj,1);
    if current_cons == 0
        CV = zeros(N,1);
    else
        CV = PopCon;
        CV = CV(:,priority);
        CV = CV(:,1:current_cons);
        CV = sum(max(0,CV),2);
    end
    
    %% Detect the dominance relation between each two solutions
    Dominate = false(N);  
    for i = 1 : N-1       
        for j = i+1 : N           
             if constraint_handing==0
                z=[PopObj CV];
                k = any(z(i,:)<z(j,:)) - any(z(i,:)>z(j,:));
                if k == 1
                   Dominate(i,j) = true;
                elseif k == -1
                   Dominate(j,i) = true;
                end
            else
                if CV(i) < CV(j)
                    Dominate(i,j) = true;
                elseif CV(i) > CV(j)
                    Dominate(j,i) = true;
                else
                    k = any(PopObj(i,:)<PopObj(j,:)) - any(PopObj(i,:)>PopObj(j,:));
                    if k == 1
                        Dominate(i,j) = true;
                    elseif k == -1
                        Dominate(j,i) = true;
                    end
                end   
            end
        end
    end
    
    %% Calculate S(i)
    S = sum(Dominate,2);
    
    %% Calculate R(i)
    R = zeros(1,N);
    for i = 1 : N
        R(i) = sum(S(Dominate(:,i)));
    end
    
    %% Calculate D(i)
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Distance = sort(Distance,2);
    D = 1./(Distance(:,floor(sqrt(N)))+2);
    
    %% Calculate the fitnesses
    Fitness = R + D';
end