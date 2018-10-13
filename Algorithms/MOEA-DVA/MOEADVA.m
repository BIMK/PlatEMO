function MOEADVA(Global)
% <algorithm> <H-N>
% A Multiobjective Evolutionary Algorithm based on Decision Variable
% Analyses for Multi-objective Optimization Problems with Large Scale
% Variables
% NCA --- 20 --- The number of sampling solutions in control variable analysis
% NIA ---  6 --- The maximum number of tries required to judge the interaction
% operator   --- DE

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [NCA,NIA] = Global.ParameterSet(20,6);
    
    %% Control variable analysis
    [DiverIndexes,ConverIndexes] = ControlVariableAnalysis(Global,NCA);
    
    %% Dividing distance variables based on two variable analyses
    [W,Global.N] = UniformPoint(Global.N,Global.M);
    [Subcomponents,Population] = DividingDistanceVariables(Global,NIA,DiverIndexes,ConverIndexes);

    %% Calculate the neighbours of each individual
    PopDec = Population.decs;
    Dis    = pdist2(PopDec(:,DiverIndexes),PopDec(:,DiverIndexes));
    Dis(logical(eye(length(Dis)))) = inf;
    [~,Neighbour] = sort(Dis,2);
    Neighbour     = Neighbour(:,1:ceil(Global.N/10));

    %% Subcomponent optimization
    if Global.M == 2; threshold = 0.01; else threshold = 0.03; end
    utility = inf;
    OldObj  = Population.objs;
    while Global.NotTermination(Population) && ~isempty(Subcomponents) && utility > threshold
        for i = 1 : length(Subcomponents)
            drawnow();
            Population = SubcomponentOptimizer(Global,Population,Neighbour,Subcomponents{i});
        end
        if ~mod(Global.gen,2)
            % For better performance on UF problems, this part should be
            % removed
            utility = sum(sum(OldObj-Population.objs))/Global.N;
            OldObj  = Population.objs;
        end
    end

    %% Uniformity optimization
    EvolveByMOEAD(Global,Population,W);
end