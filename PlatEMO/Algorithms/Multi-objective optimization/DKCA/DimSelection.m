function [dim] = DimSelection(Problem,Mask,Dec,Dec2,D,T)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yingwei Li

    Base  = zeros(1,D);
    test1 = [Base;Dec.*Mask];
    test2 = [Base;Dec2.*Mask];
    Ind   = Problem.Evaluation(test1);
    Ind2  = Problem.Evaluation(test2);
    FrontNo      = NDSort(Ind.objs,inf);
    temp_FrontNo = NDSort(Ind2.objs,inf);

    dim      = find(FrontNo<=FrontNo(1,1))-1;
    temp_dim = find(temp_FrontNo<=temp_FrontNo(1,1))-1;
    NS_num   = size(dim,2);

    if NS_num <= T  
        dim = union(temp_dim,dim);
    else
        dim = intersect(temp_dim,dim);
    end
end
