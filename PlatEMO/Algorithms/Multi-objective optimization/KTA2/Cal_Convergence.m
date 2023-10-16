function flag = Cal_Convergence(PopObj1,PopObj2,Zmin)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Zhenshou Song

    N1 = size(PopObj1,1);
    N2 = size(PopObj2,1);
    if N1~=N2
        flag = 0;
    else
        PopObj = [PopObj1;PopObj2] - repmat(Zmin,size(PopObj1,1)+size(PopObj2,1),1);
        PopObj = (PopObj)./repmat(max(PopObj,[],1) - Zmin,size(PopObj,1),1);
        Distance1 = zeros(1,N1);
        Distance2 = zeros(1,N2);
        % calculate the distance sets of CCA and CDA
        for i = 1:N1
            Distance1(i) = sqrt(sum(PopObj(i,:),2));
        end
        for i = 1: N2
            Distance2(i) = sqrt(sum(PopObj(N1+i,:),2));
        end
        % rank-sum test, alpha = 0.05
        [~,flag,~,r1,r2]=signrank_new(Distance1, Distance2,'alpha',0.05);
        if flag == 1 && (r1-r2 <0)
            flag = 0;
        end
    end
end