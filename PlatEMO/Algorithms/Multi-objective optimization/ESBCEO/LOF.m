function lof = LOF(dist)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    K=6;
    m=size(dist,1);
    distance = zeros(m,m);
    num = zeros(m,m);
    kdistance = zeros(m,1);
    count  = zeros(m,1);
    reachdist = zeros(m,m);
    lrd = zeros(m,1);
    lof = zeros(m,1);

    for i=1:m 
        [distance(i,:),num(i,:)]=sort(dist(i,:),'ascend');
        kdistance(i)=distance(i,K+1); 
        count(i) = -1;
        for j = 1:m
            if dist(i,j)<=kdistance(i)
                count(i) = count(i)+1;
            end
        end
    end

    for i = 1:m
        for j=1:i-1
            reachdist(i,j) = max(dist(i,j),kdistance(j));
            reachdist(j,i) = reachdist(i,j);
        end
    end

    for i = 1:m
        sum_reachdist=0;
        for j=1:count(i)
            sum_reachdist=sum_reachdist+reachdist(i,num(j+1));
        end
        lrd(i)=count(i)/sum_reachdist;
    end

    for i=1:m
        sumlrd=0;
        for j=1:count(i)
            sumlrd=sumlrd+lrd(num(j+1))/lrd(i);
        end
        lof(i)=sumlrd/count(i);
    end
end