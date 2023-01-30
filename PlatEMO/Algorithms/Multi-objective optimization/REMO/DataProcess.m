function [TrainIn,TrainOut,TestIn,TestOut] = DataProcess(Input,Output)
% Divide the data into the train data and test data in proportion

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    pha     = 3/4;
    index0  = find(Output==0);
    indexp1 = find(Output == 1);
    indexn1 = find(Output == -1);

    K0  = false(1,length(index0));
    Kp1 = false(1,length(indexp1));
    Kn1 = false(1,length(indexn1));

    K0(randperm(length(index0),ceil(pha*length(index0))))    = true;
    Kp1(randperm(length(indexp1),ceil(pha*length(indexp1)))) = true;
    Kn1(randperm(length(indexn1),ceil(pha*length(indexn1)))) = true;

    K        = [index0(K0);indexp1(Kp1);indexn1(Kn1)];
    TrainIn  = Input(K,:);
    TrainOut = Output(K);

    TestIn  = Input(setdiff(1:size(Input,1),K),:);
    TestOut = Output(setdiff(1:size(Input,1),K));

    Train_randindex = randperm(size(TrainOut,1),size(TrainOut,1));
    TrainIn         = TrainIn(Train_randindex,:);
    TrainOut        = TrainOut(Train_randindex);

    Test_randindex = randperm(size(TestOut,1),size(TestOut,1));
    TestIn         = TestIn(Test_randindex,:);
    TestOut        = TestOut(Test_randindex);
end