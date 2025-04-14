function C = NewSelect(P,DB,alpha,Problem)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhiyao Zhang (email: z.zhang0@csu.edu.cn)

    
    %% Preparing Data
    C     = [];
    index = [];
    for i = 1 : size(P.decs,1)
        dist2 = pdist2(real(P.decs(i,:)),real(DB.decs));
        if min(dist2) > 1e-50
            index = [index,i];
        end
    end
    if length(index) <= alpha
       PopNew = P.decs(index,:);
       if ~isempty(PopNew)
           PopNew = Problem.Evaluation(PopNew);
           C      = [C,PopNew];
       end
       return; 
    end
    
    PopDec = P.decs(index,:);
    PopObj = P.objs(index,:);
    ObjMSE = P.objmse(index,:);

    A2Obj  = DB.objs;
    [F_,~] = NDSort(A2Obj,1);
    A2Obj  = unique(A2Obj(F_==1,:),'rows');
    
    zmin   = min([A2Obj;PopObj],[],1); zmax = max([A2Obj;PopObj],[],1);
    A2Obj  = (A2Obj - zmin)./max(zmax - zmin,10e-10);
    PopObj = (PopObj - zmin)./max(zmax - zmin,10e-10);
    ObjMSE = ObjMSE./(max(zmax - zmin,10e-10).^2);
    
    %% Selection1
    [FrontNo,~] = NDSort_DIPD(PopDec,PopObj,ObjMSE,1);
    PopDec      = PopDec(FrontNo==1,:);
    PopObj      = PopObj(FrontNo==1,:);
    ObjMSE      = ObjMSE(FrontNo==1,:);
    
    if size(PopDec,1) <= alpha
        C = [C,Problem.Evaluation(PopDec)];
        return;
    end
    
    %% Selection2
    Pindex = true(1,size(PopObj,1));
    while length(find(Pindex==0)) < alpha 
        Last     = find(Pindex==1);
        Dis      = Distance(PopObj(Last,:),A2Obj);
        [~,Rank] = sort(Dis,'descend');
        PopNew   = PopDec(Last(Rank(1)),:);
        C        = [C,Problem.Evaluation(PopNew)];
        A2Obj    = [DB.objs;C.objs]; 
        [F_P,~]  = NDSort(A2Obj,1);
        A2Obj    = unique(A2Obj(F_P==1,:),'rows');
        A2Obj    = (A2Obj - zmin)./max(zmax - zmin,10e-10);
        Pindex(Last(Rank(1))) = 0;
    end
end

function dis = Distance(PopObj,OffObj)
    %% Calculate the angle-based distance between each two solutions
    dis = acos(1-pdist2(PopObj,OffObj,'cosine'));
    
    %% Calculate D
    dis = sort(dis,2);
    dis = dis(:,1);
end