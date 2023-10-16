function [XXs,Ls] = GetRelationPairs(Input,Catalog)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    C1_index = Catalog == 1;
    C2_index = Catalog ~= 1;
    C1C1 = combvec(Input(Catalog ==1,:)',Input(Catalog ==1,:)')';
    C1C2 = combvec(Input(Catalog ==1,:)',Input(Catalog ~=1,:)')';
    C2C1 = combvec(Input(Catalog ~=1,:)',Input(Catalog ==1,:)')';
    C2C2 = combvec(Input(Catalog ~=1,:)',Input(Catalog ~=1,:)')';

    t_ind     = combvec(1:sum(C1_index),1:sum(C1_index));
    t_equ_ind = t_ind(1,:) == t_ind(2,:);
    C1C1(t_equ_ind,:) = [];
    
    t_ind     = combvec(1:sum(C2_index),1:sum(C2_index));
    t_equ_ind = t_ind(1,:) == t_ind(2,:);
    C2C2(t_equ_ind,:) = [];

    t_num = ceil(size(C1C2,1)/2);
    if size(C1C1,1) > t_num && size(C2C2,1) > t_num
        C1C1 = C1C1(randperm(size(C1C1,1),t_num),:);
        C2C2 = C2C2(randperm(size(C2C2,1),t_num),:);
    elseif size(C1C1,1) < t_num
        C2C2 = C2C2(randperm(size(C2C2,1),t_num*2-size(C1C1,1)),:);
    elseif size(C2C2,1) < t_num
        C1C1 = C1C1(randperm(size(C1C1,1),t_num*2-size(C2C2,1)),:);
    end

    XXs = [C1C1;C2C2;C1C2;C2C1];
    Ls  = [zeros(size(C1C1,1),1);zeros(size(C2C2,1),1);ones(size(C1C2,1),1);-1.*ones(size(C2C1,1),1)];
end