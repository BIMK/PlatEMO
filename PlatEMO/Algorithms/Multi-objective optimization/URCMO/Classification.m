function [ flag,ll] = Classification( Population1,Population2, beita)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Kangjia Qiao

    ll = -1;

    pop1 = Population1.decs;
    conv1 = Population1.cons;
    conv1(conv1<=0)=0;
    conv1 = sum(conv1,2);
    obj1 = Population1.objs;
    pop2 = Population2.decs;
    conv2 = Population2.cons;
    conv2(conv2<=0)=0;
    conv2 = sum(conv2,2);
    obj2 = Population2.objs;

    [FrontNo,MaxFNo] = NDSort(obj1,conv1,inf);
    x1 = find(FrontNo==1);
    [FrontNo,MaxFNo] = NDSort(obj2,inf);
    x2 = find(FrontNo==1);


    if length(find(conv2(x2)<=0))==0 % When all solutions of population2 are infeasible, type-IV. Herein, 3 indicates S3
        flag = 3;
    elseif length(find(conv2(x2)>0))==0 % When all solutions of population2 are feasible, type-I.
        flag = 1;
    elseif  length(find(conv2>0))>0 && length(find(conv2>0))< size(pop1,1)  % When population2 has both feasible and infeasible solutions

        obj1 = obj1(x1,:);
        obj2 = obj2(x2,:);
        [FrontNo,MaxFNo] = NDSort([obj1;obj2],inf); 

        ll = length((find(FrontNo(1:length(x1))==1)))/length(FrontNo(1:length(x1)));  
        %  ll indicates the ratio that pop1(x1) belong to the first level in
        %  the combination of pop1(x1) and pop2(x2)

        if ll > beita
            flag = 1;
        elseif ll < 1 - beita
            flag = 3;
        else
            flag = 2;
        end
    end
end