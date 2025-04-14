function delta = Relief(Problem, Population, FrontNo)
% Relief feature selection

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    labels = zeros(1,Problem.N);
    labels(FrontNo == 1) = 1;
    features = Population.decs;
    delta    = zeros(size(features));
    distance = pdist2(features,features) + diag(inf+zeros(1,size(features,1)));
    data     = Normalize(Problem,features);
    if all(labels==1)
        delta = zeros(1,Problem.D);
    else
        for i = 1 : size(features,1)
            [~,~,nearhit_index,nearmiss_index] = SearchNear(distance,labels,i,features);
            delta(i,:) = Relevant_feature(nearhit_index,nearmiss_index,data,i);
        end
        delta = sum(delta);
    end
end

function delta = Relevant_feature(nearhit_index,nearmiss_index,data,number)
    diff_hit  = abs(data(nearhit_index,:)-data(number,:));
    diff_miss = abs(data(nearmiss_index,:)-data(number,:));
    delta     = diff_miss.^(2) - diff_hit.^(2);
end

function data = Normalize(Problem, features)
    for i = 1 : size(features,1)
        for j = 1 : size(features,2)
            data(i,j) = (features(i,j)-Problem.lower(:,j))/(Problem.upper(:,j)-Problem.lower(:,j));
        end
    end
end

function [nearhit,nearmiss,nearhit_index,nearmiss_index] = SearchNear(sample_distance,labels,number,features)
    nearhit_list  = [];
    nearmiss_list = [];
    hit_index     = [];
    miss_index    = [];
    for i = 1 : size(features, 1)
        if labels(i) == labels(number)
            nearhit_list = [nearhit_list, sample_distance(i, number)];
            hit_index    = [hit_index, i];
        else
            nearmiss_list = [nearmiss_list, sample_distance(i, number)];
            miss_index    = [miss_index, i];
        end
    end
    [~,nearhit_dis_index]  = min(nearhit_list);
    nearhit_index          = hit_index(nearhit_dis_index);
    [~,nearmiss_dis_index] = min(nearmiss_list);
    nearmiss_index         = miss_index(nearmiss_dis_index);
    nearhit                = features(nearhit_index,:);
    nearmiss               = features(nearmiss_index,:);
end