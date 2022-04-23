function [AD,Rank] = UpdateDiversityArchive(AD,Population,ND,R,Z,Xre,sigma_niche,Problem)
% Update the Diversity Archive

%------------------------------- Copyright --------------------------------
% Copyright 2017-2018 Yiping Liu
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------

    %% Non-dominated sorting
    Population = [Population,AD];
    [FrontNo,MaxFNo] = NDSort(Population.objs,ND);

    %% Select the solutions in the last front
    if length(find(FrontNo<=MaxFNo))== ND
       Next = FrontNo<=MaxFNo;
    else
       Next = FrontNo < MaxFNo ;
       NQ     = ND-sum(Next);
       Last   = find(FrontNo==MaxFNo);
       Choose = LastSelection(Population(Last).objs,Population(Last).decs,NQ,R,Z,Xre,sigma_niche,Problem); %Algorithm 4, lines 5-30
       Next(Last(Choose)) = true;
    end
    
    %% Population for next generation
    AD = Population(Next);
    Rank = FrontNo(Next);
end

function Choose = LastSelection(PopObj,PopDec,NQ,R,Z,Xre,sigma_niche,Problem)
% Select solutions with good diversity in the last front    
    N      = size(PopObj,1);
    NR     = size(R,1);
    
    %% Normalize
    PopObj = PopObj - repmat(Z,N,1);
    PopDec = (PopDec - repmat(Problem.lower,N,1))./repmat(Problem.upper-Problem.lower,N,1);
  
    %% Calculate distance between every solution and reference point in the objective space
    theta = pdist2(R,PopObj,'cosine');
    
    %% Calculate distance between every two solutions in the reminder decsion subspace
    d = pdist2(PopDec(:,Xre),PopDec(:,Xre),'chebychev');   
    
    %% Cluster
    [thmin,label] = min(theta,[],1);   
    C1 = false(NR,N);
    C2 = false(NR,N);
    for j = 1:NR
        member = find(label==j);
        member1 = label==j;
        [~,temp] = sort(thmin(member));
        for i = temp
           if any(d(member(i),and(C1(j,:),member1))<sigma_niche) 
               C2(j,member(i))=true;
           else
               C1(j,member(i))=true;
           end
        end        
    end    
    
    %% Make selected solution == NQ
    while sum(sum(C1))>NQ
        cmax = max(sum(C1,2));
        jmax = sum(C1,2)== cmax;
        temp1 = find(sum(C1(jmax,:),1)>0);        
        [~,xmax] = max(thmin(temp1));
        xmax=temp1(xmax);
        C1(label(xmax),xmax) = false;
    end   
    while sum(sum(C1))<NQ
        c2 = sum(C2,2)>0;
        cmin = min(sum(C1(c2,:),2));
        jmin = and(sum(C1,2)==cmin,c2);
        temp1 = find(sum(C2(jmin,:),1)>0);
        [~,xmin] = min(thmin(temp1));
        xmin=temp1(xmin);
        C2(label(xmin),xmin) = false;
        C1(label(xmin),xmin) = true;
    end    
    Choose = sum(C1)>0;  
end