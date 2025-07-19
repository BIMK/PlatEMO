function C = Candi_Select(PopDec,PopObj,PopCon,ObjMSE,ConMSE,Database,mu)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhiyao Zhang (email: zhiyao_zhang0@163.com)

    index = ismember(PopDec,Database.decs,'rows');
    if sum(index) == size(PopDec,1)
        C_ = [];
    else
        index_ = Selection(PopObj(~index,:),ObjMSE(~index,:),PopCon(~index,:),ConMSE(~index,:),Database,mu);
        index  = find(~index);
        C_     = PopDec(index(index_),:);
    end
    C = [];
    for i = 1 : size(C_,1)
        dist2 = pdist2(real(C_(i,:)),real(Database.decs));
        if min(dist2) > 1e-5
            C = [C;C_(i,:)];
        end
    end
end

function index = EucDistance_Select(PopObj,ALL_Obj)
    N1 = size(PopObj,1);
    N2 = size(ALL_Obj,1);
    Distance = zeros(N1,N2);

    for i = 1 : N1
        for j = 1 : N2
            Distance(i,j) = norm(PopObj(i,:)-ALL_Obj(j,:),2);    
        end
    end
    
    Distance  = sort(Distance,2);
    Distance  = Distance(:,1);
    [~,index] = max(Distance);
end

function Next = Selection(PopObj,ObjMSE,PopCon,ConMSE,Database,mu)
    %% Preparing Data
    ALL_Obj = Database.objs;
    ALL_Con = Database.cons;
    zmin    = min([ALL_Obj;PopObj]);zmax = max([ALL_Obj;PopObj]);
    ALL_Obj = (ALL_Obj - zmin )./max(zmax - zmin,10e-10);
    PopObj  = (PopObj - zmin)./max(zmax - zmin,10e-10);
    ObjMSE  =  ObjMSE./(max(zmax - zmin,10e-10).^2);

    %% Reference Set 
    global phase NI
    if phase == 2
        num = 0;
        for i= 1 : length(Database)
            if all(Database(i).cons<=0)
                num = num + 1;
            end
        end
        [FrontNo,~] = NDSort(ALL_Obj,ALL_Con,inf);
        if num > NI
            ALL_Obj = ALL_Obj(FrontNo==1,:);
        else
            i = 1;
            Next = FrontNo == i;
            while sum(Next) <= NI
                Next(FrontNo == i) = true;
                i = i + 1;
            end
            ALL_Obj = ALL_Obj(Next,:);
        end
    else
        [FrontNo,~] = NDSort(ALL_Obj,inf);
        ALL_Obj     = ALL_Obj(FrontNo==1,:);
    end

    %% Select mu points
    if phase == 2
        [FrontNo,MaxFNo] = NDSort_PDPD(PopObj,ObjMSE,PopCon,ConMSE,mu);
    else
        [FrontNo,MaxFNo] = NDSort_PDPD(PopObj,ObjMSE,mu);
    end
    Next = FrontNo < MaxFNo;
    Last = find(FrontNo == MaxFNo);
    if length(Last) == mu - sum(Next)
        Next(Last) = true;
    elseif length(Last) > mu - sum(Next)
        ALL_Obj = [ALL_Obj;PopObj(Next,:)];
        for i = 1 : mu - sum(Next)
            index             = EucDistance_Select(PopObj(Last,:),ALL_Obj);
            Next(Last(index)) = true;
            ALL_Obj           = [ALL_Obj;PopObj(Last(index),:)];
            Last(index)       = [];
        end
    end
end