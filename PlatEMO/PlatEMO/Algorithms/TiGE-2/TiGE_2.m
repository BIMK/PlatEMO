function TiGE_2(Global)
% <algorithm> <T>
% Tri-Goal Evolution Framework for CMaOPs

%------------------------------- Reference --------------------------------
% Y. Zhou, Z. Min, J. Wang, Z. Zhang, and J.Zhang, Tri-goal evolution
% framework for constrained many-objective optimization, IEEE Transactions
% on Systems Man and Cybernetics Systems, 2018.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [Epsilon0,row] = Global.ParameterSet(0.05,1.01);
    %% Generate random population
    Population = Global.Initialization();
    [fpr,fcd] = Estimation(Population.objs,1/Global.N^(1/Global.M));
    fcv = Calculate_fcv(Population); 
    Epsilon = Epsilon0;
    PopObj_1 = [fpr,fcd]; 
    [fm,~] = NDSort(PopObj_1,Global.N);
    PopObj = [fm' + Epsilon * fcv,fcv];
    [frank,~] = NDSort(PopObj,Global.N);
    fitness = frank' + fcv./(fcv+1);
    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = TournamentSelection(2,Global.N,fitness);
        Offspring  = GA(Population(MatingPool));
        [fpr,fcd] = Estimation(Offspring.objs,1/Global.N^(1/Global.M));
        fcv = Calculate_fcv(Offspring); 
        OffObj_1 = [fpr,fcd]; 
        [fm,~] =NDSort(OffObj_1,Global.N);
        OffObj = [fm' + Epsilon * fcv,fcv];
        [Population,fitness] = EnvironmentalSelection([Population,Offspring],PopObj,OffObj,Global.N);
        [fpr,fcd] = Estimation(Population.objs,1/Global.N^(1/Global.M));
        fcv = Calculate_fcv(Population);
        PopObj_1 = [fpr,fcd]; 
        [fm,~] =NDSort(PopObj_1,Global.N);
        PopObj = [fm' + Epsilon * fcv,fcv];
        Epsilon = row * Epsilon;
    end
end