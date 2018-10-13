function RPEA(Global)
% <algorithm> <O-Z>
% Many-objective Evolutionary Optimization Based on Reference Points
% alpha --- 0.4 --- Ratio of individuals being used to generate reference points
% delta --- 0.1 --- Parameter determining the difference between the reference point and the individuals

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [alpha,delta] = Global.ParameterSet(0.4,0.1);

    %% Generate random population
    Population = Global.Initialization();
    R          = GenerateRefPoints(Population,delta*(max(Population.objs,[],1)-min(Population.objs,[],1)),alpha,Global.N);
    
    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = TournamentSelection(2,Global.N,min(TchebychevDistance(Population.objs,R),[],2));
        Offspring  = Global.Variation(Population(MatingPool));
     	R          = GenerateRefPoints([Population,Offspring],delta*(max(Population.objs,[],1)-min(Population.objs,[],1)),alpha,Global.N);
        Population = EnvironmentalSelection([Population,Offspring],R,Global.N);
    end
end