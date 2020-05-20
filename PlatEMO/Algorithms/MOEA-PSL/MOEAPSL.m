function MOEAPSL(Global)
% <algorithm> <M>
% Multi-objective evolutionary algorithm based on Pareto optimal subspace
% learning

%------------------------------- Reference --------------------------------
% Y. Tian, C. Lu, X. Zhang, K. C. Tan, and Y. Jin, Solving large-scale
% multi-objective optimization problems with sparse optimal solutions via
% unsupervised neural networks, IEEE Transactions on Cybernetics, 2020.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Population initialization
    REAL = ~strcmp(Global.encoding,'binary');
    if REAL
        P   = lhsamp(Global.N,Global.D);
        Dec = P.*repmat(Global.upper-Global.lower,Global.N,1) + repmat(Global.lower,Global.N,1);
    else
        Dec = ones(Global.N,Global.D);
    end
    Mask = lhsamp(Global.N,Global.D) > 0.5;
    Population = INDIVIDUAL(Dec.*Mask);
    [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Dec,Mask,Global.N,0,0);
    
    %% Optimization
    rho = 0.5;
    while Global.NotTermination(Population)
        Site = rho > rand(1,ceil(Global.N/2));
        if any(Site)
            [rbm,dae,allZero,allOne] = ModelTraining(Mask(FrontNo==1,:),Dec(FrontNo==1,:),REAL);
        else
            [rbm,dae,allZero,allOne] = deal([],[],[],[]);
        end
        MatingPool = TournamentSelection(2,ceil(Global.N/2)*2,FrontNo,-CrowdDis);
        [OffDec,OffMask] = Operator(Dec(MatingPool,:),Mask(MatingPool,:),rbm,dae,Site,allZero,allOne,Global.lower,Global.upper,REAL);
        Offspring = INDIVIDUAL(OffDec.*OffMask);
        [Population,Dec,Mask,FrontNo,CrowdDis,sRatio] = EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Global.N,length(Population),2*sum(Site));
        rho = (rho+sRatio)/2;
    end
end