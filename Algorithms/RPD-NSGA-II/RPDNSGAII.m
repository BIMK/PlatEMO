function RPDNSGAII(Global)
% <algorithm> <R>
% Reference point dominance-based NSGA-II

%------------------------------- Reference --------------------------------
% M. Elarbi, S. Bechikh, A. Gupta, L. B. Said, and Y. S. Ong, A new
% decomposition-based NSGA-II for many-objective optimization, IEEE
% Transactions on Systems, Man, and Cybernetics: Systems, 2018, 48(7):
% 1191-1210.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate the reference points and random population
    [RPSet,Global.N] = UniformPoint(Global.N,Global.M);
    Population       = Global.Initialization();
    [~,FrontNo,d2]   = EnvironmentalSelection(Population,RPSet,Global.N);

    %% Optimization
    while Global.NotTermination(Population) 
        MatingPool = TournamentSelection(2,Global.N,FrontNo,d2);
        Offspring  = GA(Population(MatingPool));
        [Population,FrontNo,d2] = EnvironmentalSelection([Population,Offspring],RPSet,Global.N);
    end
end