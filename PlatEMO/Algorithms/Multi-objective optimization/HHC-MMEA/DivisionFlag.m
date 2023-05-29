function [subp1,subp2,DivisionFlag] = DivisionFlag(Mask,FrontNo)
% Calculate the farthest two solutions and their similarity

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    mask = Mask(FrontNo==1,:);
    dis = pdist2(mask,mask,'hamming');
    M = max(max(dis));%maximum distance
    [leader1,leader2] = find(dis==M,1);
    [~,subleader1]=ismember(mask(leader1,:),Mask,'rows');
    [~,subleader2]=ismember(mask(leader2,:),Mask,'rows');
    subp1 = subleader1;
    subp2 = subleader2;
    s = sum(mask(leader1,:)&mask(leader2,:))/min(sum(mask(leader1,:)),sum(mask(leader2,:)));%相似度
    if M>0.05 && s<=0.5
        DivisionFlag = 1;
        %% hierarchical clustering(divide operation)
        for i = 1:size(Mask,1)
            if i~=subleader1 && i~=subleader2
                dis1 = pdist2(mask(leader1,:),Mask(i,:),'hamming');
                dis2 =  pdist2(mask(leader2,:),Mask(i,:),'hamming');
                if dis1 <= dis2
                    subp1 = [subp1,i];
                else
                    subp2 = [subp2,i];
                end
            end
        end
    else
        DivisionFlag = 0;
    end
    if size(Mask,1) <50
        DivisionFlag = 0;
    end
end