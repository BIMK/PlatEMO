function SPEAR(Global)
% <algorithm> <O-Z>
% A Strength Pareto Evolutionary Algorithm Based on Reference Direction for
% Multi-objective and Many-objective Optimization

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate the reference directions (general approach)
    [W,Global.N] = UniformPoint(Global.N,Global.M);
    % Largest acute angle between two neighbouring reference directions
    cosine = 1 - pdist2(W,W,'cosine');
    cosine(logical(eye(length(cosine)))) = 0;
    theta  = max(min(acos(cosine),[],2));
    
    %% Generate random population
    Population = Global.Initialization();

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = MatingSelection(Population.objs,20);
        Offspring  = Global.Variation(Population([1:Global.N,MatingPool]),Global.N);
        QObj       = ObjectiveNormalization([Population,Offspring]);
        [Ei,Angle] = Associate(QObj,W);
        FV         = FitnessAssignment(Ei,QObj,Angle,theta);
        Population = EnvironmentalSelection(Global,[Population,Offspring],Ei,FV);
    end
end