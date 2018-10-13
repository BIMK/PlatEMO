function Two_Arch2(Global)
% <algorithm> <O-Z>
% An Improved Two-Archive Algorithm for Many-Objective Optimization
% CAsize --- --- Convergence archive size
% p      --- --- The parameter of fractional distance
% operator   --- EAreal

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
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
        Offspring = [Global.Variation(ParentC,inf,@EAreal,{[],[],0,0}),...
                     Global.Variation(ParentM,inf,@EAreal,{0,0,[],[]})];
        CA = UpdateCA(CA,Offspring,CAsize);
        DA = UpdateDA(DA,Offspring,Global.N,p);
    end
end