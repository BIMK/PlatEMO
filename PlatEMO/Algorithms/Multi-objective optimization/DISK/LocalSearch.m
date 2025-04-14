function DB = LocalSearch(P,W,ideal,wmax,Model,DB,Problem)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhiyao Zhang (email: z.zhang0@csu.edu.cn)

    w     = 1;
    [N,~] = size(P.decs);
    while w <= wmax   
        OffDec1 = OperatorGA(Problem,P.decs);
        OffDec2 = OperatorDE_current_rand_1(Problem,P.decs);
        OffDec3 = OperatorDE_rand_1(Problem,P.decs);
        OffDec4 = OperatorDE_current_rand_1(Problem,P.decs);
        P.decs  = [P.decs;OffDec1;OffDec2;OffDec3;OffDec4];
        P.decs  = unique(P.decs,'rows');
        
        [P.objs,P.objmse] = model_predict(Model,P.decs);
        P.objmse = sqrt(P.objmse);
                
        fitness = max(abs(P.objs - ideal).*W,[],2) - 2*mean(P.objmse,2);       
        
        [~,Rank] = sort(fitness);    
        P.decs   = P.decs(Rank(1:N),:);
        P.objs   = P.objs(Rank(1:N),:);
        P.objmse = P.objmse(Rank(1:N),:);
        w = w + 1;
    end
    
    fitness  = max(abs(P.objs - ideal).*W,[],2) - 2*mean(P.objmse,2);
    [~,Rank] = sort(fitness);
    PopNew   = P.decs(Rank(1),:); 
    dist2    = pdist2(real(PopNew),real(DB.decs));
    if min(dist2) > 1e-50
        DB = [DB,Problem.Evaluation(PopNew)];
    end
end

function [OffObj,Off_ObjMSE] = model_predict(Model,OffDec)
    N          = size(OffDec,1);
    Len_obj    = length(Model);
    OffObj     = zeros(N,Len_obj);
    Off_ObjMSE = zeros(N,Len_obj);
   
    for i = 1 : N
        for j = 1 : Len_obj
            [OffObj(i,j),~,Off_ObjMSE(i,j)] = predictor(OffDec(i,:),Model{j});
        end
    end
    OffObj     = real(OffObj);
    Off_ObjMSE = abs(real(Off_ObjMSE));
end