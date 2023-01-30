function [Fitness, extreme] = CalFitness(PopObj)
% Calculate the fitness of each solution
% F(x) = f(x) + d(x) Algorithm 3

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Li Li

    [N,M]  = size(PopObj); %N  No. of solutions, M No. of Objectives
   
   %% q
    %Zmax = max(PopObj,[],1); % Identify the max point
    Zmin = min(PopObj,[],1); % Identify the ideal point
    % Identify the extreme points
    W = zeros(M) + 1e-6;
    W(logical(eye(M))) = 1;
    ASF = zeros(N,M);
    for i = 1 : M
        ASF(:,i) = max((PopObj-repmat(Zmin,N,1))./repmat(W(i,:),N,1),[],2);
    end
    [~,extreme] = min(ASF,[],1); %'extreme' is the extreme solutions
    % Calculate the intercepts
    Hyperplane = PopObj(extreme,:)\ones(M,1); % linear equation X=A\B???A*X=B???hyperplane is the solution of linear equation PopObj(extreme,:)x Hyperplane = ones(M,1)
    a = (1./Hyperplane)';
    if any(isnan(a))
        a = max(PopObj,[],1);
    end
    % Normalization
    PopObj = (PopObj-repmat(Zmin,N,1))./repmat(a-Zmin,N,1);
    
    %---------------w_nad
    w_nad = a; % set a as the nadir point

    Zmin = min(PopObj,[],1); % Identify the ideal point
    for i = 1 : M
        ASF_q(:,i) = max((PopObj-repmat(Zmin,N,1))./repmat(w_nad,N,1),[],2);
    end
    [~,key_point] = min(ASF_q,[],1);
    %q
    q = sqrt(sum(PopObj(key_point).^2,2)) * sqrt(M); 
    if q <= 1
        fx = sum(PopObj-repmat(Zmin,N,1),2);
    else
        fx = max(PopObj-repmat(Zmin,N,1),[],2);
    end
  
    %% Calculate the DMD between each two solutions  
    Distance = inf(N);
    for i = 1 : N
        for j = [1:i-1,i+1:N]
            denominator = 1 + min(PopObj(i,:),PopObj(j,:)); 
            numerator = max(PopObj(i,:),PopObj(j,:)) - min(PopObj(i,:),PopObj(j,:));  
            dis = numerator./ denominator;
            Distance(i,j) = sum(dis);
        end
    end
      
    %% Calculate D(i)
    Distance = sort(Distance,2);
    k = floor(sqrt(N)); %k-NN
    d = 1./(Distance(:,k)+2); % dimensionality margin distance
    
    %% Calculate the fitnesses
    Fitness = fx+d;
end