function U = Operator(lu,X,bestX)
% Differential evolution

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

     V = mutation(X,bestX,0.8,2);
     U = crossover(lu,X,V,1,1);
end

function V = mutation(X,bestX,F,mutationStrategy)
    NP = size(X,1);
    V  = zeros(size(X));
    for i = 1 : NP
        nrandI = 5;
        r = randi([1,NP],1,nrandI);
        for j = 1 : nrandI
            equalr(j) = sum(r==r(j));
        end
        equali   = sum(r==i);
        equalval = sum(equalr) + equali;
        while equalval > nrandI
            r = randi([1,NP],1,nrandI);
            for j = 1 : nrandI
                equalr(j) = sum(r==r(j));
            end
            equali   = sum(r==i);
            equalval = sum(equalr)+equali;
        end
        switch mutationStrategy
            case 1
                %mutationStrategy=1:DE/rand/1;
                V(i,:) = X(r(1),:) + F*(X(r(2),:)-X(r(3),:));
            case 2
                %mutationStrategy=2:DE/best/1;
                V(i,:) = bestX + F*(X(r(1),:)-X(r(2),:));
            case 3
                %mutationStrategy=3:DE/rand-to-best/1;
                V(i,:) = X(i,:) + F*(bestX-X(i,:)) + F*(X(r(1),:)-X(r(2),:));
            case 4
                %mutationStrategy=4:DE/best/2;
                V(i,:) = bestX + F*(X(r(1),:)-X(r(2),:)) + F*(X(r(3),:)-X(r(4),:));
            case 5
                %mutationStrategy=5:DE/rand/2;
                V(i,:) = X(r(1),:) + F*(X(r(2),:)-X(r(3),:)) + F*(X(r(4),:)-X(r(5),:));
        end
    end
end

function U = crossover(lu,X,V,CR,crossStrategy)
    [NP,D] = size(X);
    U      = zeros(NP,D);
    switch crossStrategy
        case 1
            for i = 1 : NP
                jRand = floor(rand*D);
                for j = 1 : D
                    if rand<CR || j==jRand
                        U(i,j) = V(i,j);
                    else
                        U(i,j) = X(i,j);
                    end
                end
            end
        case 2
            for i = 1 : NP
                j = floor(rand*D);
                L = 0;
                U = X;
                U(i,j) = V(i,j);
                j = mod(j+1,D);
                L = L + 1;
                while rand<CR && L<D
                    U(i,j) = V(i,j);
                    j = mod(j+1,D);
                    L = L + 1;
                end
            end
    end
    for i = 1 : NP
        for j = 1 : D
            minin = min(lu(1,j),lu(2,j));
            maxax = max(lu(1,j),lu(2,j));
            while U(i,j)>maxax || U(i,j)<minin
                U(i,j) = minin + rand*(maxax-minin);
            end
        end
    
    end
end