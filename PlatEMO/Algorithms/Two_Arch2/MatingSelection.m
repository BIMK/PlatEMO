function [ParentC,ParentM] = MatingSelection(CA,DA,N)
% The mating selection of Two_Arch2

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    CAParent1 = randi(length(CA),1,ceil(N/2));
    CAParent2 = randi(length(CA),1,ceil(N/2));
    Dominate  = any(CA(CAParent1).objs<CA(CAParent2).objs,2) - any(CA(CAParent1).objs>CA(CAParent2).objs,2);  
    ParentC   = [CA([CAParent1(Dominate==1),CAParent2(Dominate~=1)]),...
                 DA(randi(length(DA),1,ceil(N/2)))];
    ParentM   = CA(randi(length(CA),1,N));
end