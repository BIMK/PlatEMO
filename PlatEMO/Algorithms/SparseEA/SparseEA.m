function SparseEA(Global)
% <algorithm> <S>
% Evolutionary algorithm for sparse multi-objective optimization problems

%------------------------------- Reference --------------------------------
% Y. Tian, X. Zhang, C. Wang, and Y. Jin, An evolutionary algorithm for
% large-scale sparse multi-objective optimization problems, IEEE
% Transactions on Evolutionary Computation, 2019.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Population initialization
    % Calculate the fitness of each decision variable
    TDec    = [];
    TMask   = [];
    TempPop = [];
    Fitness = zeros(1,Global.D);
    REAL    = ~strcmp(Global.encoding,'binary');
    for i = 1 : 1+4*REAL
        if REAL
            Dec = unifrnd(repmat(Global.lower,Global.D,1),repmat(Global.upper,Global.D,1));
        else
            Dec = ones(Global.D,Global.D);
        end
        Mask       = eye(Global.D);
        Population = INDIVIDUAL(Dec.*Mask);
        TDec       = [TDec;Dec];
        TMask      = [TMask;Mask];
        TempPop    = [TempPop,Population];
        Fitness    = Fitness + NDSort([Population.objs,Population.cons],inf);
    end
    % Generate initial population
    if REAL
        Dec = unifrnd(repmat(Global.lower,Global.N,1),repmat(Global.upper,Global.N,1));
    else
        Dec = ones(Global.N,Global.D);
    end
    Mask = zeros(Global.N,Global.D);
    for i = 1 : Global.N
        Mask(i,TournamentSelection(2,ceil(rand*Global.D),Fitness)) = 1;
    end
    Population = INDIVIDUAL(Dec.*Mask);
    [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection([Population,TempPop],[Dec;TDec],[Mask;TMask],Global.N);

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool       = TournamentSelection(2,2*Global.N,FrontNo,-CrowdDis);
        [OffDec,OffMask] = Operator(Dec(MatingPool,:),Mask(MatingPool,:),Fitness,REAL);
        Offspring        = INDIVIDUAL(OffDec.*OffMask);
        [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Global.N);
    end
end