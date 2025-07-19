function C = Candi_Select(PopDec,PopObj,PopCon,ObjMSE,ConMSE,Database,mu,Problem)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhiyao Zhang (email: zhiyao_zhang0@163.com)

    %% Preparing Data
    index = ismember(PopDec,Database.decs,'rows');
    if sum(index) == size(PopDec,1)
        C = [];
        return;
    elseif sum(~index) <= mu
        C_ = PopDec(~index,:);
        C  = [];
        for i = 1 : size(C_,1)
            dist2 = pdist2(real(C_(i,:)),real(Database.decs));
            if min(dist2) > 1e-5
                C = [C;C_(i,:)];
            end
        end
        return;
    end
    
    PopDec = PopDec(~index,:);
    PopObj = PopObj(~index,:); ObjMSE  = ObjMSE(~index,:);
    PopCon = PopCon(~index,:); ConMSE  = ConMSE(~index,:);

    A2Obj  = Database.objs;
    A2Con  = Database.cons;
    zmin   = min([A2Obj;PopObj],[],1); zmax = max([A2Obj;PopObj],[],1);
    A2Obj  = (A2Obj - zmin )./max(zmax - zmin,10e-10);
    PopObj = (PopObj - zmin)./max(zmax - zmin,10e-10);
    ObjMSE = ObjMSE./(max(zmax - zmin,10e-10).^2);
    
    PopCon = max(PopCon,0);A2Con = max(A2Con,0);
    PopCon = PopCon./max([PopCon;A2Con],[],1);
    A2Con  = A2Con./max([PopCon;A2Con],[],1);
    PopCV  = sum(PopCon,2);A2CV = sum(A2Con,2);
    PopCV  = PopCV./max(max([PopCV;A2CV],10e-10));
    A2CV   = A2CV./max(max([PopCV;A2CV],10e-10));
    
    %% Reference Set 
    num = length(find(all(A2Con<=0,2)));
    [FrontNo,~] = NDSort(A2Obj,A2Con,inf);
    if num >= Problem.N
        A2Obj = A2Obj(FrontNo==1,:);
        A2CV  = A2CV(FrontNo==1,:);
    else
        i    = 1;
        Next = FrontNo == i;
        while sum(Next) < Problem.N
            Next(FrontNo == i) = true;
            i = i + 1;
        end
        A2Obj = A2Obj(Next,:);
        A2CV  = A2CV(Next,:);
    end

    %% Select mu points
    [FrontNo,MaxFNo] = NDSort_CPPD(PopObj,ObjMSE,PopCon,ConMSE,mu);
    Next = FrontNo < MaxFNo;
    Last = find(FrontNo == MaxFNo);

    if length(Last) == mu - sum(Next)
        Next(Last) = true;
    elseif length(Last) > mu - sum(Next)
        A2Obj = [A2Obj;PopObj(Next,:)];
        A2CV  = [A2CV;PopCV(Next,:)];
        for i = 1 : mu - sum(Next)
            index = EucDistance_Select([PopObj(Last,:),PopCV(Last,:)],[A2Obj,A2CV]);
            Next(Last(index)) = true;
            A2Obj = [A2Obj;PopObj(Last(index),:)];
            A2CV  = [A2CV;PopCV(Last(index),:)];
            Last(index)=[];
        end
    end
    
    %% Output
    C_ = PopDec(Next,:);
    C  = [];
    for i = 1 : size(C_,1)
        dist2 = pdist2(real(C_(i,:)),real(Database.decs));
        if min(dist2) > 1e-5
            C = [C;C_(i,:)];
        end
    end
end

function index = EucDistance_Select(PopObj,ALL_Obj)
    N1       = size(PopObj,1);
    N2       = size(ALL_Obj,1);
    Distance = zeros(N1,N2);

    %% Calculate the shifted distance between each two solutions
    for i = 1 : N1
        SPopObj = max(ALL_Obj,repmat(PopObj(i,:),N2,1));
        for j = 1 : N2
            Distance(i,j) = norm(PopObj(i,:)-SPopObj(j,:),2);    
        end
    end
    Distance  = sort(Distance,2);
    [~,index] = max(Distance(:,1));
end