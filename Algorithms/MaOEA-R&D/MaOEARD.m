function MaOEARD(Global)
% <algorithm> <H-N>
% Many-Objective Evolutionary Algorithm: Objective Space Reduction +
% Diversity

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate the weight vectors
    W = zeros(Global.M) + 1e-6;
    W(logical(eye(Global.M))) = 1;
    
    %% Generate random population
    Global.N   = max(ceil(Global.N/Global.M)*Global.M,2*Global.M);
    Population = Global.Initialization();
    
    %% Objective space reduction
    % Reducing the objective space
    while Global.NotTermination(Population) && Global.evaluated < Global.evaluation/2
        % Classification
        [Subpopulation,Z] = Classification(Population,W);
        % Evolve each subpopulation
        for i = 1 : Global.M
            MatingPool = randi(Global.N/Global.M,1,Global.N/Global.M);
            Offspring  = Global.Variation(Subpopulation{i}(MatingPool));
            Subpopulation{i} = [Subpopulation{i},Offspring];
            ASF = max((Subpopulation{i}.objs-repmat(Z,length(Subpopulation{i}),1))./repmat(W(i,:),length(Subpopulation{i}),1),[],2);
            [~,rank] = sort(ASF);
            Subpopulation{i} = Subpopulation{i}(rank(1:Global.N/Global.M));
        end
        Population = [Subpopulation{:}];
    end
    % Identify the target points
    TP = UpdateTP(Population,W);

    %% Diversity improvement
    % Generate random population around TPs
    PopDec     = repmat(TP.decs,Global.N/Global.M-1,1);
    PopDec     = PopDec.*(1+randn(size(PopDec))/5);
    Population = [TP,INDIVIDUAL(PopDec)];
    % Evolve
    while Global.NotTermination(Population)
        MatingPool = randi(Global.N,1,Global.N);
        Offspring  = Global.Variation(Population(MatingPool));
        Population = EnvironmentalSelection([Population,Offspring],Global.N,TP.objs);
        TP         = UpdateTP([Population,TP],W);
    end
end