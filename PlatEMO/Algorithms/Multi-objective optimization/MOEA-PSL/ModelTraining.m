function [rbm,dae,allZero,allOne] = ModelTraining(Mask,Dec,REAL)
% Training RBM and DAE

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Determine the size of hidden layers
    allZero = all(~Mask,1);
    allOne  = all(Mask,1);
    other   = ~allZero & ~allOne;
    K       = sum(mean(abs(Mask(:,other).*Dec(:,other))>1e-6,1)>rand(1,sum(other)));
    K       = min(max(K,1),size(Mask,1));
    
    %% Train RBM and DAE
    rbm = RBM(sum(other),K,10,1,0,0.5,0.1);
    rbm.train(Mask(:,other));
    if REAL
        dae = DAE(size(Dec,2),K,10,size(Dec,1),0.5,0.5,0.1);
        dae.train(Dec);
    else
        dae = [];
    end
end