function [SourceId,TF1] = SourceTaskrand(Problem,SubPopulation,Dec,Mask,FrontNo,CrowdDis,EachN,Fitness,FitnessDec,pDecIndex,Buquan)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    TaskNum  = size(Problem.SubM,2);
    SourceId = zeros(1,TaskNum);
    TF1 = [];
    for i = 1 : TaskNum
        TF1 = [TF1;Fitness{i}/5,zeros(1,Buquan(i))];
        TF1(TF1==0) = max(Fitness{i}/5)+1;
    end
    mse_matrix = zeros(TaskNum);  
    for i = 1 : TaskNum
        for j = 1:TaskNum
            mse_matrix(i, j) = mean((TF1(i, :) - TF1(j, :)).^2); 
        end
        arr          = find(mse_matrix(i,:) ~= 0);
        shuffled_arr = arr(randperm(length(arr)));
        SourceId(i)  = shuffled_arr(1);
    end        
end