function MOEADVA(Global)
% <algorithm> <M>
% Multi-objective evolutionary algorithm based on decision variable
% analyses
% NCA --- 20 --- The number of sampling solutions in control variable analysis
% NIA ---  6 --- The maximum number of tries required to judge the interaction

%------------------------------- Reference --------------------------------
% X. Ma, F. Liu, Y. Qi, X. Wang, L. Li, L. Jiao, M. Yin, and M. Gong, A
% multiobjective evolutionary algorithm based on decision variable analyses
% for multiobjective optimization problems with large-scale variables, IEEE
% Transactions Evolutionary Computation, 2016, 20(2): 275-298.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
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
    while Global.NotTermination(Population)
        for i = 1 : length(Subcomponents)
            drawnow();
            Population = SubcomponentOptimizer(Population,Neighbour,Subcomponents{i});
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