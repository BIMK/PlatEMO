function Population = EnvironmentalSelection(OffSpring,W,N)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming

    S       = []; 	% S is the set used for output (replaced by population as output)
    Sc      = [];  	% Sc is used to collect feasible solutions
    PopObj  = OffSpring.objs;
    PopCon  = OffSpring.cons;
    Fitness = CalFitness1(PopObj,PopCon);
    Sc=[Sc,OffSpring(Fitness<1)];

    if length(Sc)==N
        Population=Sc;
    elseif length(Sc)>N
        [FrontNO,MaxNO]=NDSort(Sc.objs,inf);
        if MaxNO==1
           %% ¦È-dominated sorting
            FrontNO = tNDSort(Sc.objs,W);
            MaxNO=max(FrontNO);
            for i=1:MaxNO
                S=cat(2,S,Sc(FrontNO==i));
                if length(S)>=N
                    break;
                end
            end
        else
            for i=1:MaxNO
                S=cat(2,S,Sc(FrontNO==i));
                if length(S)>=N
                    break;
                end
            end
        end
        while length(S)>N
            %normalization
            Zmax=max(S.objs,[],1);
            Zmin=min(S.objs,[],1);
            SPopObj=(S.objs-repmat(Zmin,size(S.objs,1),1))./(repmat(Zmax,size(S.objs,1),1)-repmat(Zmin,size(S.objs,1),1));
            
            [~,Region] = max(1-pdist2(SPopObj,W,'cosine'),[],2);% associate each solution in S with their corresponding subregion
            [value,~]=sort(Region,'ascend');
            flag=max(value);
            counter=histc(value,1:flag);                        % counter denotes the number of indiviudals in each subregion
            [~,most_crowded]=max(counter);
            S_crowdest=S(Region==most_crowded);                 % S_crowdest is the set of individuals from the most crowded subregion
            dist=pdist2(S_crowdest.objs,S_crowdest.objs);
            dist(dist==0)=inf;
            [row,~]=find(min(min(dist))==dist);
            St=S_crowdest(row);                                 % St is the set of individuals having the smallest distance in S_crowdest
            [~,Region_St] = max(1-pdist2(St.objs,W,'cosine'),[],2);
            Z = min(St.objs,[],1);
            g_tch=max(abs(St.objs-repmat(Z,length(St),1))./W(Region_St,:),[],2);
            [~,order]=max(g_tch);
            x_wrost=St(order);
            S=setdiff(S,x_wrost);
        end
        Population=S;
    elseif length(Sc)<N
        SI=setdiff(OffSpring,Sc);                                      % SI is the set of infeasible solutions in Hc
        f1=sum(max(0,SI.cons),2);
        [~,Region_SI] = max(1-pdist2(SI.objs,W,'cosine'),[],2);
        Z = min(SI.objs,[],1) ;
        f2=max(abs(SI.objs-repmat(Z,length(SI),1))./W(Region_SI,:),[],2);
        PopObj=[f1,f2];
        [FrontNO,MaxNO]=NDSort(PopObj,inf);
        S=[S,Sc];
        for i=1:MaxNO
            S=cat(2,S,SI(FrontNO==i));
            if length(S)>=N
                last=i;
                break;
            end
        end
        F_last=SI(FrontNO==last);                               % find the individuals in the last front joined into S
        delete_n=size(S,2)-N;
        CV=sum(max(0,F_last.cons),2);
        [~,index]=sort(CV,'descend');
        x_wrost=F_last(index(1:delete_n));
        S=setdiff(S,x_wrost);
        Population=S;
    end
end