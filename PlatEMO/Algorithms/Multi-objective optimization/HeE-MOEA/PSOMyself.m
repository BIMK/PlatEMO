function index = PSOMyself(x, y, D)  
% Select input features for the generation of different inputs by PSO

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    m=20;% Population Size
    iter=30;% Number of Generation
    threshold=0.5;c1=2;c2=2;w=1.05;% Parameter of CSO
    
    % Limit of Position
    lu = [zeros(1, D); ones(1, D)];
    XRRmin = repmat(lu(1, :), m, 1);
    XRRmax = repmat(lu(2, :), m, 1);
    VelMax=(XRRmax-XRRmin)/10;

    % Position
    p = rand(m, D);
    p=XRRmin + (XRRmax - XRRmin) .*p;
    p(find(p<threshold))=0;p(find(p>=threshold))=1;
    v = VelMax .* rand(m,D);

    fitness = PSO_CostFunction(p, x, y);

    [bestever, index] = min(fitness);
    zbest=p(index,:);fitnesszbest=bestever;
    gbest=p;fitnessgbest=fitness;

    for i = 1 : iter

        v = w*v + c1*rand(m,D).*(gbest-p) + c2*rand(m,D).*(repmat(zbest, m, 1)-p);
        v = max(min(v,VelMax), -VelMax);
        p = p + v;
        p(find(p<threshold))=0;p(find(p>=threshold))=1;

        fitness = PSO_CostFunction(p, x, y); 

        index = find(fitness < fitnessgbest);
        fitnessgbest(index) = fitness(index);
        gbest(index,:) = p(index,:);
        [bestever, index] = min(fitness);
        if  bestever < fitnesszbest
            fitnesszbest = bestever;
            zbest = p(index, :);
        end
    end
    index=find(zbest==1);
end

function fitness=PSO_CostFunction(p, x, y)
    m=size(p,1);
    fitness=100*ones(m,1);
    for i=1:m
        index=find(p(i,:)==1);
        if ~isempty(index)
            xx=x(:,index);
            D=Discor(xx, y);
            R=Discor(xx,[]);
            fitness(i)=-0.8*D+0.2*R;
        end
    end
end

function D=Discor(x, y)
    [n, x_n]=size(x);
    if ~isempty(y)
        a=norm(x);
        b=norm(y);
        A=a-repmat(mean(a,2),1,n)-repmat(mean(a,1),n,1)+mean(mean(a));
        B=b-repmat(mean(b,2),1,n)-repmat(mean(b,1),n,1)+mean(mean(b));
        covxy=sum(sum(A.*B))/(n^2);
        covx=sum(sum(A.^2))/(n^2);
        covy=sum(sum(B.^2))/(n^2);
        D=sqrt(covxy)/sqrt(sqrt(covx)*sqrt(covy));

    else
        if x_n==1
            D=0;
        else
            for i=1:x_n
                xa=x(:,i);
                index=setdiff([1:x_n],i);
                xb=x(:,index);    
                a=norm(xa);
                b=norm(xb);
                A=a-repmat(mean(a,2),1,n)-repmat(mean(a,1),n,1)+mean(mean(a));
                B=b-repmat(mean(b,2),1,n)-repmat(mean(b,1),n,1)+mean(mean(b));
                covxy=sum(sum(A.*B))/(n^2);
                covx=sum(sum(A.^2))/(n^2);
                covy=sum(sum(B.^2))/(n^2);
                D(1,i)=sqrt(covxy)/sqrt(sqrt(covx)*sqrt(covy));
            end
            D=mean(D);
        end
    end
end