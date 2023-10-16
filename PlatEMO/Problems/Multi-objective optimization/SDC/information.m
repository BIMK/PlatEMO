function a = information(index)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Kangjia Qiao (email: qiaokangjia@yeah.net)

    CEC_problem = [1,2,3,6,9,10,11,12,14,18,19,24,15,5,1];
    shape_problem = [1,2,1,2,1,2,1,2,1,2,1,2,2,1,1];
    b = [10,100,15,115,19,125,10,100,15,115,19,115,125,15,10];
    Distance_problem = [2,1,4,4,3,5,5,3,3,2,1,3,5,1,2]; 

    HCT = 0.5;
    high_type= [1,2,1,1,2,2,2,2,1,1,1,2,2,2,1];
    DCT = 0.5;
    dis_type = [2,1,1,2,2,2,1,1,2,2,1,1,1,1,2];

    i = index;
    a = [CEC_problem(i),Distance_problem(i),...
        HCT,high_type(i),DCT,dis_type(i),shape_problem(i),b(i)];
end