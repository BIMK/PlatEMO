function UpdatedDA=UpdateDA(CA,DA,Q,W)
% Update DA
% CA is the Archive that has  been updated.
% DA is the Archive that has not been updated
% Q  is the set of offspring
% W  is the set of weight vectors

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    S=[];                                                       % S is the set used for output
    Hd=[DA,Q];
    N=size(W,1);
    [~,Region_Hd] = max(1-pdist2(real(Hd.objs),W,'cosine'),[],2);
    [~,Region_CA] = max(1-pdist2(CA.objs,W,'cosine'),[],2);     % associat the individuals in Hd and CA with their corresponding subregions
    itr=1;
    while length(S)<N
        for i=1:N                                               % here i denotes the order of the subregion
            current_c=find(Region_CA==i);                       % current_c denotes that the current_c_th individual is/are already in the ith region
            if length(current_c)<itr 
                for j=1:itr-length(current_c)                   % j denotes the number of solutions from Hd that need to join into the region(i)
                    current_d=find(Region_Hd==i);
                    if isempty(current_d)~=1
                        [FrontNO,~]=NDSort(Hd(current_d).objs,inf);
                        O=Hd(current_d(FrontNO==1));            % O is the set of nondominated solutions from region(i) in Hd
                        [~,Region_O] = max(1-pdist2(real(O.objs),W,'cosine'),[],2);
                        Z = min(O.objs,[],1);
                        g_tch=max(abs(O.objs-repmat(Z,length(O),1))./W(Region_O,:),[],2);
                        [~,order]=min(g_tch);
                        x_best=O(order); 

                        Hd(current_d(Hd(current_d)==x_best))=[];% update Region_Hd
                        if isempty(Hd)==1
                            Region_Hd=[];
                        else
                            [~,Region_Hd] = max(1-pdist2(real(Hd.objs),W,'cosine'),[],2);
                        end    
                        if length(S)<N                          % add the best individual into S
                            S=cat(2,S,x_best);
                        end
                    else
                        break;
                    end
               end
            end
            if length(S)==N
                break;
            end
        end
        itr=itr+1;
    end
    UpdatedDA=S;
end