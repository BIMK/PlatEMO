function LMEA(Global)
% <algorithm> <H-N>
% A Decision Variable Clustering Based Evolutionary Algorithm for
% Large-scale Many-objective Optimization
% nSel ---  5 --- Number of selected solutions for decision variable clustering
% nPer --- 50 --- Number of perturbations on each solution for decision variable clustering
% nCor ---  5 --- Number of selected solutions for decision variable interaction analysis

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [nSel,nPer,nCor] = Global.ParameterSet(5,50,5);

    %% Generate random population
    Population = Global.Initialization();

    %% Detect the group of each distance variable
    [PV,DV] = VariableClustering(Global,Population,nSel,nPer);
    DVSet   = CorrelationAnalysis(Global,Population,DV,nCor);

    %% Optimization
    while Global.NotTermination(Population)
        % Convergence optimization
        for i = 1 : 10
            drawnow();
        	Population = ConvergenceOptimization(Global,Population,DVSet);
        end
        % Distribution optimization
        for i = 1 : Global.M
            drawnow();
            Population = DistributionOptimization(Global,Population,PV);
        end
    end
end

function Population = ConvergenceOptimization(Global,Population,DVSet)
    N   = length(Population);
    Con = calCon(Population.objs);
    for i = 1 : length(DVSet)
        for j = 1 : length(DVSet{i})
            % Select parents
            MatingPool = TournamentSelection(2,2*N,Con);
            % Generate offsprings
            OffDec = Population.decs;
            if isequal(Global.operator,@DE)
                NewDec = Global.VariationDec(Population([1:N,MatingPool]).decs,N,@DE,{[],[],Global.D/length(DVSet{i})/2,[]});
            else
                NewDec = Global.VariationDec(Population(MatingPool).decs,N,@EAreal,{[],[],Global.D/length(DVSet{i})/2,[]});
            end
            OffDec(:,DVSet{i}) = NewDec(:,DVSet{i});
            Offspring          = INDIVIDUAL(OffDec);
            % Update each solution
            allCon  = calCon([Population.objs;Offspring.objs]);
            Con     = allCon(1:N);
            newCon  = allCon(N+1:end);
            updated = Con > newCon;
            Population(updated) = Offspring(updated);
            Con(updated)        = newCon(updated);
        end
    end
end

function Population = DistributionOptimization(Global,Population,PV)
% Distribution optimization

    N            = length(Population);
    OffDec       = Population(TournamentSelection(2,N,calCon(Population.objs))).decs;
    NewDec       = Global.VariationDec(Population(randi(N,1,N)).decs,N,@EAreal);
    OffDec(:,PV) = NewDec(:,PV);
    Offspring    = INDIVIDUAL(OffDec);
    Population   = EnvironmentalSelection([Population,Offspring],N);
end

function Con = calCon(PopuObj)
% Calculate the convergence of each solution

    FrontNo = NDSort(PopuObj,inf);
    Con     = sum(PopuObj,2);
    Con     = FrontNo'*(max(Con)-min(Con)) + Con;
end