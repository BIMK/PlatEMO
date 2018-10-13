function [TrainIn,TrainOut,TestIn,TestOut] = DataProcess(Input,Output)
% Divide the data into the train data and test data in proportion

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He

        index1 = find(Output>0.5);
        index0 = find(Output<=0.5);
        K1     = false(1,length(index1));
        K0     = false(1,length(index0));
        K1(randperm(length(index1),ceil(3/4*length(index1)))) = true;
        K0(randperm(length(index0),ceil(3/4*length(index0)))) = true;
        K        = [index1(K1);index0(K0)];
        TrainIn  = Input(K,:);
        TrainOut = Output(K);
        TestIn   = Input(setdiff(1:size(Input,1),K),:);
        TestOut  = Output(setdiff(1:size(Input,1),K));
end