function rNSGAII(Global)
% <algorithm> <O-Z>
% The r-Dominance: A New Dominance Relation for Interactive Evolutionary
% Multicriteria Decision Making
% Points ---     --- Set of preferred points
% W      ---     --- Set of weight vector for each preferred point
% delta  --- 0.1 --- Non-r-dominance threshold

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [Points,W,delta] = Global.ParameterSet(zeros(1,Global.M)+0.5,ones(1,Global.M),0.1);

    %% Generate random population
    Population = Global.Initialization();
    FrontNo    = NrDSort(Population.objs,inf,Points,W,1-(1-delta)*Global.gen/Global.maxgen);
    CrowdDis   = CrowdingDistance(Population.objs,FrontNo);

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
        Offspring  = Global.Variation(Population(MatingPool));
        [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Global.N,Points,W,1-(1-delta)*Global.gen/Global.maxgen);
    end
end