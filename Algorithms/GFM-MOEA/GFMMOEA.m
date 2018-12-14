function GFMMOEA(Global)
% <algorithm> <G>
% Generic front modeling based multi-objective evolutionary algorithm
% theta --- 0.2 --- Penalty parameter
% fPFE  --- 0.1 --- Frequency of employing generic front modeling

%------------------------------- Reference --------------------------------
% Y. Tian, X. Zhang, R. Cheng, C. He, and Y. Jin, Guiding evolutionary
% multi-objective optimization with generic front modeling, IEEE
% Transactions on Cybernetics, 2018.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [theta,fPFE] = Global.ParameterSet(0.2,0.1);
    
    %% Generate random population
    Population = Global.Initialization();
    FrontNo    = NDSort(Population.objs,inf);
    zmin       = min(Population.objs,[],1);
    % Calculate the fitness of each solution
    [P,A]     = deal(ones(1,Global.M));
    [App,Dis] = CalFitness(Population.objs-repmat(zmin,length(Population),1),P,A);
    Dis       = sort(Dis,2);
    Crowd     = Dis(:,1) + 0.1*Dis(:,2);

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = TournamentSelection(2,Global.N,FrontNo,-theta*App-(1-theta)*Crowd);
        Offspring  = GA(Population(MatingPool));
        zmin       = min([zmin;Offspring.objs],[],1);
        if ~mod(Global.gen,ceil(fPFE*Global.maxgen)) || fPFE == 0
            [P,A] = GFM(Population(FrontNo==1).objs-repmat(zmin,sum(FrontNo==1),1));
        end
        [Population,FrontNo,App,Crowd] = EnvironmentalSelection([Population,Offspring],P,A,zmin,theta,Global.N);
    end
end