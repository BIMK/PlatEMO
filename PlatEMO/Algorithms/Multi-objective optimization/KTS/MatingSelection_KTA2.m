function [ParentCobj,ParentCdec,ParentMobj,ParentMdec] = MatingSelection_KTA2(CAobj,CAdec,DAobj,DAdec,N)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhenshou Song

    CAParent1  = randi(size(CAobj,1),1,ceil(N/2));
    CAParent2  = randi(size(CAobj,1),1,ceil(N/2));
    Dominate   = any(CAobj(CAParent1,:)<CAobj(CAParent2,:),2) - any(CAobj(CAParent1,:)>CAobj(CAParent2,:),2);  
    ParentCobj = [CAobj([CAParent1(Dominate==1),CAParent2(Dominate~=1)],:);...
                 DAobj(randi(size(DAobj,1),1,ceil(N/2)),:)];
    ParentCdec = [CAdec([CAParent1(Dominate==1),CAParent2(Dominate~=1)],:);...
                 DAdec(randi(size(DAdec,1),1,ceil(N/2)),:)];
    ParentMobj = CAobj(randi(size(CAobj,1),1,N),:);
    ParentMdec = CAdec(randi(size(CAdec,1),1,N),:);
end