function [Population,FrontNo,crd] = ESelection(Population,N,Ruq,Rnode,theta,Zmin)
% The environmental selection of DEA-GNG

%--------------------------------------------------------------------------
% Copyright Yiping Liu
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------

    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,N);
    Next = FrontNo < MaxFNo;
    F1 = Population(FrontNo==1);
    Zmax = max(F1.objs,[],1); % Nadir Point
    
    %% Select the solutions in the last front
    Last   = find(FrontNo==MaxFNo);
    [Choose,crd] = LastSelection(Population(Next).objs,Population(Last).objs,N-sum(Next),Ruq,Rnode,theta,Zmin,Zmax);
    Next(Last(Choose)) = true;
    
    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);        
end

function [Choose, crd] = LastSelection(PopObj1,PopObj2,K,Ruq,Rnode,theta,Zmin,Zmax)
% Select part of the solutions in the last front

    %% Initialize
    PopObj = [PopObj1;PopObj2]; %Candidate solutions
    R = [Ruq;Rnode]; % Reference Vectors   
    [N,~]  = size(PopObj);
    N1     = size(PopObj1,1);
    N2     = size(PopObj2,1);
    NR     = size(R,1);
    NR1    = size(Ruq,1);   
        
    %% Normalization [0-1]
    PopObj = (PopObj - repmat(Zmin,N,1))./repmat(Zmax-Zmin,N,1);
    
    %% Scalarizing function value
    g = zeros(1,N); 
    d1 = zeros(1,N); %d1 in PBI
    theta = [Inf(1,NR1),theta]; %theta in PBI
        
    %% Associate each solution with one reference vector
    % Calculate the distance of each solution to each reference vector
    Cosine   = 1 - pdist2(PopObj,R,'cosine');
    NormP = sqrt(sum(PopObj.^2,2));
    Distance2 = repmat(NormP,1,NR).*sqrt(1-Cosine.^2);    
    % Associate each solution with its nearest reference vector
    [d2,pi] = min(Distance2',[],1);  
    
    %% Calculate the number of associated solutions except for the last front of each reference point
    rho = hist(pi(1:N1),1:NR);
       
    %% Calculate scalarizing functions (PBI)
    for i = N1+1:N
        d1(i) = NormP(i).*Cosine(i,pi(i));
        if theta(pi(i))==Inf
            g(i) = d2(i);
        else         
            g(i) = d1(i)+theta(pi(i)).*d2(i); 
        end    
    end
       
    %% Select solutions
    Choose  = false(1,N2);
    Zchoose = true(1,NR);
    while sum(Choose) < K
        % Select the least crowded reference point
        Temp = find(Zchoose);      
        Jmin = find(rho(Temp)==min(rho(Temp)));
        j = Temp(Jmin(randi(length(Jmin)))); 
        I = find(Choose==0 & pi(N1+1:end)==j);       
        if ~isempty(I)
            if rho(j) == 0
                [~,s] = min(g(N1+I));
            else
                s = randi(length(I));
            end
            Choose(I(s)) = true;
            rho(j) = rho(j) + 1;
        else
            Zchoose(j) = false;
        end
    end   
    crd = rho(pi);
    crd = [crd(1:N1),crd(N1+find(Choose))];   
end