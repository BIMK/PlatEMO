function [H,dist]=ShanNonx(train,data,m,D,k)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    H = zeros(size(data,1),1);
    for i = 1 : size(data,1)
        [knn,dist] = nbselect(train(:,D+1:D+m),data(i,D+1:D+m),'K',k,m);
        sum=zeros(1,D);
        Px=zeros(size(knn,2),D);
        for p=1:D
            for q=1:size(knn,2)
                sum(p)=sum(p)+train(knn(q),p);
            end
        end
        for p = 1 : D
            for q = 1 : size(knn,2)                
                if sum(p)==0
                   Px(q,p)=1/size(knn,2);
                else
                   Px(q,p)=train(knn(q),p)/sum(p);
                end
                if p==D
                    if q==size(knn,2)
                        Px=NdSort(Px,size(knn,2),D);
                    end
                end
                 if p==D
                     if Px(q,p)==0
                         h=0;
                     else
                         h=Px(q,p)*log2(Px(q,p));
                     end
                     H(i)=H(i)-h;
                 end
            end
        end
    end

end

function P=NdSort(Px,q,D)
    [~,I]=sort(Px);
    for i=1:q/2
        for j=1:D
            Px=swap(Px,I,j,i,q-i+1);
        end
    end
    P=Px;
end

function P=swap(Px,I,j,i,p)
    t=Px(I(i,j),j);
    Px(I(i,j),j)=Px(I(p,j),j);
    Px(I(p,j),j)=t;
    P=Px;
end

function [idx,dist] = nbselect(fitness,part,varargin)
    if varargin{1} == 'K'
        k = varargin{2};
        [idx,dist] = knnsearch(fitness(:,1:end),part(:,1:end),'Distance','euclidean','NSMethod','kdtree','K',k);  
    end
end