function Two_Arch2(Global)
% <algorithm> <T>
% Two-archive algorithm 2
% CAsize --- --- Convergence archive size
% p      --- --- The parameter of fractional distance

%------------------------------- Reference --------------------------------
% H. Wang, L. Jiao, and X. Yao, Two_Arch2: An improved two-archive
% algorithm for many-objective optimization, IEEE Transactions on
% Evolutionary Computation, 2015, 19(4): 524-541.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [CAsize,p] = Global.ParameterSet(Global.N,1/Global.M);
    
    %% Generate random population
    Population = Global.Initialization();
    CA = UpdateCA([],Population,CAsize);
    DA = UpdateDA([],Population,Global.N,p);
    
    %% Optimization
    while Global.NotTermination(DA)
        [ParentC,ParentM] = MatingSelection(CA,DA,Global.N);
        Offspring         = [GA(ParentC,{1,20,0,0}),GA(ParentM,{0,0,1,20})];
        CA = UpdateCA(CA,Offspring,CAsize);
        DA = UpdateDA(DA,Offspring,Global.N,p);
    end
end