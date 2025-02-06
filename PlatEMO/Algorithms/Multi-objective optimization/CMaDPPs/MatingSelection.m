function ParentC = MatingSelection(CA,DA,N,z,znad)
% The mating selection of Two_Arch2

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming (email: 20151000334@cug.edu.cn)

    EA=[CA,DA];
    CAObj=EA.objs;
    [N1,~]=size(CAObj);
    CAObj2 = (CAObj-repmat(z,N1,1))./(repmat(znad,N1,1)-repmat(z,N1,1));
    D = pdist2(CAObj2,CAObj2,'cosine');
    D=D+eye(N1);
    [cos,mincos]=min(D);
    CAO=sum(CAObj2.^2,2);
    minCAO=CAO(mincos);
    ch=(minCAO-CAO)>0;
    ch3=(1-ch).*mincos'+(ch).*(1:N1)';
    cos=1-cos;
    cos=(cos-min(cos))./repmat((max(cos)-min(cos)),1,N1);
    ParentC=[];
    for i=1:1:2*N
        k=randi(N1);
        if rand<cos(k)
            ParentC=[ParentC,EA(ch3(k))];
        else
            ParentC=[ParentC,EA(k)];
        end
    end
end