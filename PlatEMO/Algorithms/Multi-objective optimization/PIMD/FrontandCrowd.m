function [FrontNo,CrowdDis] = FrontandCrowd(objs,p)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yang Li (email: liyangwust@163.com)

    %% Round the objective values
    objs = round(objs, 6);
    
    %% Non-dominated sorting
    [nInd,~] = size(objs);
    [FrontNo,MaxFNo] = NDSort(objs,nInd);
    
    front1 = objs(FrontNo==1,:);
    [front1,normalization,IdealPoint,Extreme] = NormCalc(front1);
    
    CrowdDis = zeros(1,nInd);
    
    [m,~] = size(front1);
    F1CrowdDis = zeros(1,m) ;
    F1CrowdDis(Extreme) = Inf;
    selected = false(1,m);
    selected(Extreme) = true;
    nn = vecnorm(front1, p, 2);
    distances = pdist2(front1, front1, 'minkowski', p);
    distances = distances ./ repmat(nn, 1, m);
    neighbors = 2;
    remaining = 1:m;
    remaining = remaining(~selected);
    for i = 1 : m-sum(selected)-1
        maxim = mink(distances(remaining, selected),neighbors,2);
        [d, index] = max(sum(maxim,2));
        best = remaining(index);
        remaining(index) = [];
        selected(best)=true;
        F1CrowdDis(1,best) = d;
    end
    CrowdDis(FrontNo==1) = F1CrowdDis;
    
    for i = 2 : MaxFNo
        front = objs(FrontNo==i,:);
        [m,~] = size(front);
        front = front./repmat(normalization',m,1);
        CrowdDis(FrontNo==i) = 1./pdist2(front, IdealPoint, 'minkowski', p);
    end
end

