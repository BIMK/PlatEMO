function Choose = SelectionMSE(Obj,MSE,RealFirstObj,N,Z,Zmin)
% Nodominated sorting considering the uncertainty of the points, the
% sorting rule is K1,L1,K2,L2....

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    MSE2 = prod(MSE,2);
    MSE3 = MSE2.^(1/size(MSE,2));
    all  = [Obj,MSE3];
    if isempty(Zmin)
        Zmin = ones(1,size(Z,2));
    end

    [NRealFront,M] = size(RealFirstObj);
    NPop = size(Obj,1);
    
    R_S = zeros(NRealFront,NPop);      
    for i = 1 : NRealFront
        x        = sum(repmat(RealFirstObj(i,:),NPop,1)-Obj<=0,2) == M;  
        R_S(i,x) = 1;  
    end      
    
    index=(sum(R_S,1)==0);
    NonDomPopAll = all(index,:);
    OtherPopAll  = all(index==0,:);
   
    [FrontNo1,~] = NDSort(NonDomPopAll,N);
    OtherPopAll(:,M+1:end) = -1.*OtherPopAll(:,M+1:end);
    [FrontNo2,~] = NDSort(OtherPopAll,N);
    nonDom    = find(index==1);
    other     = find(index==0);
    Choose    = false(1,size(index,2));
    NonDomPop = Obj(index,:);
    OtherPop  = Obj(index==0,:);
    tnum      = 0;
    i = 0;
    while tnum < N
        i = i + 1;
        current1 = FrontNo1 == i;
        if tnum+sum(current1) <= N
            Choose(nonDom(current1)) = true ;
            tnum = tnum+sum(current1);
        else
            temp1 = false(1,size(FrontNo1,2));
            Last1 = find(FrontNo1==i);
            C1    = LastSelection(Obj(Choose,:),NonDomPop(current1,:),N-sum(Choose),Z,Zmin);
            temp1(Last1(C1)) = true;
            Choose(nonDom(temp1)) = true ;           
            tnum = tnum + sum(C1);
        end
        if tnum == N
            break;
        elseif tnum > N
            fprintf('error')
        else
            current2 = FrontNo2 == i;
            if tnum+sum(current2) <= N
                Choose(other(current2)) = true ;
                tnum = tnum + sum(current2);
            else
                temp2 = false(1,size(FrontNo2,2));
                Last2 = find(FrontNo2==i);
                C2    = LastSelection(Obj(Choose,:),OtherPop(current2,:),N-sum(Choose),Z,Zmin);
                temp2(Last2(C2))     = true;
                Choose(other(temp2)) = true ;           
                tnum = tnum+sum(C2);
            end    
        end
    end    
end

function Choose = LastSelection(PopObj1,PopObj2,K,Z,Zmin)
% Select part of the solutions in the last front

    PopObj = [PopObj1;PopObj2] - repmat(Zmin,size(PopObj1,1)+size(PopObj2,1),1);
    [N,M]  = size(PopObj);
    N1     = size(PopObj1,1);
    N2     = size(PopObj2,1);
    NZ     = size(Z,1);

    %% Normalization
    % Detect the extreme points
    Extreme = zeros(1,M);
    w       = zeros(M)+1e-6+eye(M);
    for i = 1 : M
        [~,Extreme(i)] = min(max(PopObj./repmat(w(i,:),N,1),[],2));
    end
    if size(unique(Extreme),1)~=M
        a = max(PopObj,[],1)';
    else
        Hyperplane = PopObj(Extreme,:)\ones(M,1);
        a = 1./Hyperplane;       
    end

    % Normalization
    PopObj = PopObj./repmat(a',N,1);
    
    %% Associate each solution with one reference point
    % Calculate the distance of each solution to each reference vector
    Cosine   = 1 - pdist2(PopObj,Z,'cosine');
    Distance = repmat(sqrt(sum(PopObj.^2,2)),1,NZ).*sqrt(1-Cosine.^2);
    % Associate each solution with its nearest reference point
    [d,pi] = min(Distance',[],1);

    %% Calculate the number of associated solutions except for the last front of each reference point
    rho = hist(pi(1:N1),1:NZ);
    
    %% Environmental selection
    Choose  = false(1,N2);
    Zchoose = true(1,NZ);
    % Select K solutions one by one
    while sum(Choose) < K
        % Select the least crowded reference point
        Temp = find(Zchoose);
        Jmin = find(rho(Temp)==min(rho(Temp)));
        j    = Temp(Jmin(randi(length(Jmin))));
        I    = find(Choose==0 & pi(N1+1:end)==j);
        % Then select one solution associated with this reference point
        if ~isempty(I)
            if rho(j) == 0
                [~,s] = min(d(N1+I));
            else
                s = randi(length(I));
            end
            Choose(I(s)) = true;
            rho(j) = rho(j) + 1;
        else
            Zchoose(j) = false;
        end
    end
end