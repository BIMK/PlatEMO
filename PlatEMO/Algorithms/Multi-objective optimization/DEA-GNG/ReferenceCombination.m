function [Ruq,net] = ReferenceCombination(Ru,net)
% Combine uniform reference vectors and nodes in GNG (Algorihtm 4)

%--------------------------------------------------------------------------
% Copyright Yiping Liu
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------

    numNode = size(net.NodeS,1);    
    net.NodeP = net.NodeS;
    
    %% Map nodes to hyperplane
    for i = 1:numNode
       x = 1./sum(net.NodeS(i,:));
       net.NodeP(i,:) = net.NodeS(i,:)*x; 
    end    
    
    %% Average distance ammong nodes
    Distance1 = pdist2(net.NodeP,net.NodeP).*net.edge;   
    AvgDis = sum(Distance1(:))./sum(net.edge(:)); 
    
    %% Min distance among Ru
    Distance2 = pdist2(Ru,Ru);
    Distance2(Distance2 == 0) = 1;   
    MinDis = min(min(Distance2));  
  
    %% Choose the smaller one
    AvgDis = min(AvgDis,MinDis);

    %% Remove some reference vectors in Ru which are too close to nodes
    Distance3 = pdist2(Ru,net.NodeP);      
    Choose = all(Distance3 > AvgDis,2)==1;   
    Ruq = Ru(Choose,:);   
end