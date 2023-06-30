function CrowdDis = Crowding(Pop)
% Harmonic average distance of each solution in the decision space
% return: the crowding distance of each individual

    [N, ~]=size(Pop);
    K=N-1;
    Z = min(Pop,[],1);
    Zmax = max(Pop,[],1);
    pop=(Pop-repmat(Z,N,1))./repmat(Zmax-Z,N,1);
    distance=pdist2(pop,pop);
    [value,~]=sort(distance,2);
    CrowdDis=K./sum(1./value(:,2:N),2);
end