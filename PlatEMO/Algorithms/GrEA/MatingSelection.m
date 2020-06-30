function MatingPool = MatingSelection(PopObj,div)
% The mating selection of GrEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N = size(PopObj,1);

    %% Calculate the grid location of each solution
    fmax = max(PopObj,[],1);
    fmin = min(PopObj,[],1);
    lb   = fmin-(fmax-fmin)/2/div;
    ub   = fmax+(fmax-fmin)/2/div;
    d    = (ub-lb)/div;
    lb   = repmat(lb,N,1);
    d    = repmat(d,N,1);
    GLoc = floor((PopObj-lb)./d); 
    GLoc(isnan(GLoc)) = 0;
    
    %% Calculate the GD value of each solution
    GD = zeros(N)+inf;
    for i = 1 : N-1
        for j = i+1 : N
            GD(i,j) = sum(abs(GLoc(i,:)-GLoc(j,:)));
            GD(j,i) = GD(i,j);
        end
    end
    
    %% Calculate the GCD value of each solution
    GD  = max(size(PopObj,2)-GD,0);
    GCD = sum(GD,2);
    
    %% Binary tournament selection
    Parents1   = randi(N,1,N);
    Parents2   = randi(N,1,N);
    Dominate   = any(PopObj(Parents1,:)<PopObj(Parents2,:),2) - any(PopObj(Parents1,:)>PopObj(Parents2,:),2);
    GDominate  = any(GLoc(Parents1,:)<GLoc(Parents2,:),2) - any(GLoc(Parents1,:)>GLoc(Parents2,:),2);
    MatingPool = [Parents1(Dominate==1 | GDominate==1),...
                  Parents2(Dominate==-1 | GDominate==-1),...
                  Parents1(Dominate==0 & GDominate==0 & GCD(Parents1)<=GCD(Parents2)),...
                  Parents2(Dominate==0 & GDominate==0 & GCD(Parents1)>GCD(Parents2))];
end