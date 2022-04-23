function theta = TunePBI(net,eps)
% Adapt scalarizing functions by tune theta in PBI

%--------------------------------------------------------------------------
% Copyright Yiping Liu
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------
 
    N = size(net.NodeS,1);
    theta=zeros(1,N);
    for i=1:N
        Neighbor = find(net.edge(i,:)==1);
        N1 = length(Neighbor);
        if N1==0
            theta(i)=Inf;
            continue;
        end
        EdgeVector = repmat(net.NodeS(i,:),N1,1)-net.NodeS(Neighbor,:);      
        flag=0;
        for j=1:N1
            if all(EdgeVector(j,:)==0)
                flag=1;
                break;
            end
        end
        if flag==1
            theta(i)=Inf;
            continue; 
        end
        Cosine = 1 - pdist2(net.NodeS(i,:),EdgeVector,'cosine');
        Cosine = max(Cosine);
        Angle = acos(Cosine);
        Angle =  Angle - eps;
        if Angle<0
            Angle=0;
        end
        theta(i)=1./tan(Angle);
        if theta(i)<0
            theta(i)=0;
        end
    end   
end