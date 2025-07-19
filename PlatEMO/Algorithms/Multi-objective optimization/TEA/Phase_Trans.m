function [phase,ct] = Phase_Trans(A2,C,ct,ct_max,phase)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhiyao Zhang (email: zhiyao_zhang0@163.com)

    if phase == 1
        feasible_num   = 0;
        feasible_index = [];
        index          = 0;
        for i= 1 : length(A2)
            if all(A2(i).cons<=0)
                feasible_num   = feasible_num + 1;
                feasible_index = [feasible_index,i];
            end
        end
        
        if feasible_num ~= 0
            feasible_newindex   = [];feasible_newnum = 0;
            infeasible_newindex = [];
            for i= 1 : length(C)
                if all(C(i).cons<=0)
                    feasible_newnum   = feasible_newnum + 1;
                    feasible_newindex = [feasible_newindex,i];
                else
                    infeasible_newindex = [infeasible_newindex,i];
                end
            end
           
            if ((isempty(feasible_newindex)==0 && set_dominate(C(feasible_newindex),A2(feasible_index)) == 3) || isempty(feasible_newindex)) &&...
               (isempty(infeasible_newindex)||(isempty(infeasible_newindex)==0 && (set_dominate(C(infeasible_newindex),A2(feasible_index)) == 3 ||...
               set_dominate(C(infeasible_newindex),A2(feasible_index)) == 1)))
                ct = ct + 1;
                if ct >= ct_max
                    index = 1;
                end
            else
                ct = 0;
            end
        end
        
        if (feasible_num ~= 0) && (phase == 1) && (index == 1)
            phase = 2;
        else
            phase = 1;
        end
    end
end

function flag = set_dominate(A,B)
    %flag = 1: one of A dominte B
    %flag = 2: B dominte A
    %flag = 3: A is nondominated with B
    %flag = 4: other

    [FrontNo,~]    = NDSort(B.objs,inf);
    B              = B(FrontNo == 1);
    [FrontNo,~]    = NDSort(A.objs,inf);
    A              = A(FrontNo == 1);
    Aobj           = A.objs;
    Bobj           = B.objs;
    Asize          = length(A);
    Bsize          = length(B);
    dominate_index = zeros(Asize,Bsize + 1);
    for i = 1 : Asize
        for j = 1 : Bsize
            if all(Aobj(i,:) == Bobj(j,:))
                dominate_index(i,j) = 3;
            elseif all(Aobj(i,:) <= Bobj(j,:))
                dominate_index(i,j) = 1;
            elseif all(Aobj(i,:) >= Bobj(j,:))
                dominate_index(i,j) = 2;
            else
                dominate_index(i,j) = 3;
            end
        end

        uni = unique(dominate_index(i,1:Bsize));
        uni = sort(uni);
        if length(uni) == 1
            dominate_index(i,Bsize + 1) = uni;
        elseif length(uni) == 2
            if all(uni == [1,3])
                dominate_index(i,Bsize + 1) = 1;
            elseif all(uni == [2,3])
                dominate_index(i,Bsize + 1) = 2;
            else
                dominate_index(i,Bsize + 1) = 4;
            end
        else
            dominate_index(i,Bsize + 1) = 4;
        end
    end

    uni = unique(dominate_index(:,Bsize + 1))';

    if length(uni) == 1
        flag = uni;
    elseif length(uni) == 2
        if all(uni == [1,3])
            flag = 1;
        elseif all(uni == [2,3])
            flag = 2;
        else
            flag = 4;
        end
    else
        flag = 4;
    end
end