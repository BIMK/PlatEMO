function Population= OCPO_Selection(OCPO,pop_num,Problem)  
%The environment selection method of OCPO

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhiqiang Zeng (email: zhiqiang.zeng@outlook.com)

    [N1,~]       = size(OCPO.objs);
    [FrontNo1,~] = NDSort(OCPO.objs,OCPO.cons,inf);
    CrowdDis1    = CrowdingDistance(OCPO.objs,FrontNo1);
    [~,r1]       = sortrows([FrontNo1',-CrowdDis1']);
    Rc(r1)       = 1 : N1;  % Ranking based on CDP


    [FrontNo2,~] = NDSort(OCPO.objs,0,inf);
    CrowdDis2    = CrowdingDistance(OCPO.objs,FrontNo2);
    [~,r2]       = sortrows([FrontNo2',-CrowdDis2']);
    Rp(r2)       = 1 : N1;  % Ranking based on non dominated sorting

    a     = 1-Problem.FE/Problem.maxFE;
    pro_l = 1-length(find(sum(max(0,OCPO.cons),2)>0))/N1;   % Calculate feasibility rate
    b     = 1.0/((Problem.FE/Problem.maxFE)^2+1)-0.5;


    r        = a*(1-b)+pro_l*b;
    R_sum    = (1-r)*Rc+r*Rp;   % To fuse two rankings.
    [~,Rank] = sort(R_sum);

     FEA = OCPO(find(FrontNo1==1));
    
     PopObj = FEA.objs;
     [P, M] = size(PopObj);
     if P > pop_num
         if Problem.M <= 2
             Del                    = Truncation(PopObj, P - pop_num);
             remained_solution      = 1 : P;
             remained_solution(Del) = [];
             FEA                    = FEA(remained_solution);
         else
             [W, N]   = UniformPoint(pop_num, M);
             zmin     = min(PopObj, [], 1) - 1e-6;
             zmax     = max(PopObj, [], 1);
             PopObj   = (PopObj - zmin) ./ (zmax - zmin);
             Distance = sqrt(sum(PopObj.^2, 2));
             Angle_Pop_to_W = sin(acos(1 - pdist2(W, PopObj, 'cosine')));
             S   = FEA;
             FEA = [];
             for i = 1 : N
                 Fitness   = Distance' .* Angle_Pop_to_W(i,:);
                 [~,index] = min(Fitness);
                 FEA       = [FEA,S(index)];
             end
             if N < pop_num
                 FEA = [FEA,OCPO(Rank(1:pop_num-N))];
             end
        end
        Population = FEA;
    else
        Population = OCPO(Rank(1:pop_num));
    end
 end

 function Del = Truncation(PopObj,K)
% Select part of the solutions by truncation

    %% Truncation
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,size(PopObj,1));
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
 end