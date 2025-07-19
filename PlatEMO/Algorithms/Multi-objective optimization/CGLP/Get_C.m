function [Mm,hisPareto] = Get_C(hisPareto,Mm,T)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    POF       = hisPareto{T}.F';
    POS       = hisPareto{T}.X';
    tempPof   = POF;
    tempPos   = POS;
    [~,index] = sort(POF(:,1));
    for i = 1 : length(index)
        POF(i,:) = tempPof(index(i),:);
        POS(i,:) = tempPos(index(i),:);
    end
    hisPareto{T}.F = POF';
    hisPareto{T}.X = POS';
    
    curr_FrontNo  = NDSort(hisPareto{T}.F',inf);
    curr_NDS1     = [hisPareto{T}.X(:,find(curr_FrontNo==1));find(curr_FrontNo==1)]';
    curr_NDF1     = [hisPareto{T}.F(:,find(curr_FrontNo==1));find(curr_FrontNo==1)]';
    curr_CrowdDis = CrowdingDistance(curr_NDF1(:,1:size(tempPof,2)));
    [~,sort_curr] = sort(curr_CrowdDis','descend');
    curr_NDS      = curr_NDS1(sort_curr,1:size(tempPos,2));
    Mm{T}         = curr_NDS;
end