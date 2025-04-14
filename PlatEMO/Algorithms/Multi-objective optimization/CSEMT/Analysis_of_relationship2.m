function flag = Analysis_of_relationship2(temp1,temp2,target1,target2,beta)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Kangjia Qiao (email: qiaokangjia@yeah.net)

    conv1 = temp1.cons;
    conv1(conv1<=0) = 0;
    obj1  = temp1.objs;

    conv2 = temp2.cons;
    conv2(conv2<=0) = 0;
    obj2 = temp2.objs;
    
    FrontNo = NDSort(obj1,conv1(:,target1),inf);
    x1      = find(FrontNo==1);
    FrontNo = NDSort(obj2,conv2(:,target2),inf);
    x2      = find(FrontNo==1);
    
    obj1    = obj1(x1,:);
    obj2    = obj2(x2,:);
    FrontNo = NDSort([obj1;obj2],inf);

    alpha1 = length(find(FrontNo(1:length(x1))==1))/length(x1);
    alpha2 = length(find(FrontNo(length(x1)+1:end)==1))/length(x2);
    
    if alpha1 >= beta && alpha2 < beta
        flag = 0;
    elseif  alpha1 < beta && alpha2 >= beta
        flag = 1;
    elseif  alpha1 < beta && alpha2 < beta
        flag = 2;
    elseif  alpha1 >= beta && alpha2 >= beta
        flag = 3;
    end
end