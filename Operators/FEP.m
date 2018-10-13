function Offspring = FEP(Global,Parent)
% <operator> <real>
% Faster evolutionary programming

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    ParentDec = Parent.decs;
    [N,D]     = size(ParentDec);
    ParentEta = Parent.adds(rand(N,D));
    
    %% FEP
    tau  = 1/sqrt(2*sqrt(D));
    tau1 = 1/sqrt(2*D);
    GaussianRand  = repmat(randn(N,1),1,D);
    GaussianRandj = randn(N,D);
    CauchyRandj   = trnd(1,N,D);
    NewDec = ParentDec + ParentEta.*CauchyRandj;
    NewEta = ParentEta.*exp(tau1*GaussianRand+tau*GaussianRandj);
    
    Offspring = INDIVIDUAL(NewDec,NewEta);
end