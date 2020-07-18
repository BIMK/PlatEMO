function DGEA(Global)
% <algorithm> <D>
% A direction guided algorithm for large-scale multiobjective optimization
% RefNo     ---  10 --- Number of reference vectors for offspring generation
% operation ---   0 --- Operation of the environmental selection
%------------------------------- Reference --------------------------------
% C. He, R. Cheng, and D. Yazdani, Adaptive offspring generation for
% evolutionary large-scale multiobjective optimization, IEEE Transactions
% on System, Man, and Cybernetics: Systems, 2020.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He

    %% Generate random population
    [operation, RefNo]= Global.ParameterSet(1, 20);
	[V,Global.N] = UniformPoint(Global.N,Global.M);
    Population   = Global.Initialization();
    Offspring    = Global.Initialization();
    Arc = Population;
    
    %% Optimization
    while Global.NotTermination(Arc)
        [Population,FrontNo] = PreSelection([Population,Offspring],V,(Global.gen/Global.maxgen)^2,RefNo);
        Offspring = DirectionReproduction(Global,Population,FrontNo,RefNo);
        switch operation
            case 0
                Arc = [Population,Offspring]; % Without selection
            case 1 
                Arc = subRVEA([Arc,Offspring],V,(Global.gen/Global.maxgen)^2);
            case 2
                Arc = subNSGAII([Arc,Offspring],Global.N);
            case 3 
                Arc = subIBEA([Arc,Offspring],Global.N,0.05);
            case 4 
                Arc = subSPEA2([Arc,Offspring],Global.N);
        end
    end
end
