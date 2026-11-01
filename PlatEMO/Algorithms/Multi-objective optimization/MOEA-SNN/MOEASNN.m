classdef MOEASNN < ALGORITHM
% <2026> <multi> <real/binary> <large/none> <sparse>
% Multi-objective evolutionary algorithm based on spiking neural networks

%------------------------------- Reference --------------------------------
% F. Ming, X. Shi, W. Pan, and Y. Jin. Spike-driven sparse binary
% representation learning for large-scale multi-objective optimization.
% IEEE Transactions on Evolutionary Computation, 2026.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming (email: mingfei@westlake.edu.cn)

    methods
        function main(Algorithm,Problem)
            %% Population initialization
            P          = UniformPoint(Problem.N,Problem.D,'Latin');
            Dec        = P.*repmat(Problem.upper-Problem.lower,Problem.N,1) + repmat(Problem.lower,Problem.N,1);
            Dec(:,Problem.encoding==4) = 1;
            P_Dec      = InitialPop(Problem);
            Mask       = P_Dec.decs;
            Population = Problem.Evaluation(Dec.*Mask);
            [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Dec,Mask,Problem.N,0,0);

            %% Optimization
            rho = 0.5;
            while Algorithm.NotTerminated(Population)
                Site = rho > rand(1,ceil(Problem.N/2));
                if any(Site)
                    [snn,dae,allZero,allOne] = ModelTraining(Mask,Dec,any(Problem.encoding~=4));
                else
                    [snn,dae,allZero,allOne] = deal([]);
                end
                MatingPool       = TournamentSelection(2,ceil(Problem.N/2)*2,FrontNo,-CrowdDis);
                [OffDec,OffMask] = Operator(Problem,Dec(MatingPool,:),Mask(MatingPool,:),snn,dae,Site,allZero,allOne);
                Offspring        = Problem.Evaluation(OffDec.*OffMask);
                [Population,Dec,Mask,FrontNo,CrowdDis,sRatio] = EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Problem.N,length(Population),2*sum(Site));
                rho              = (rho+sRatio)/2;
            end
        end
    end
end

function Pop = InitialPop(Problem)
    N = Problem.N;
    D = Problem.D;
    T = min(D, N * 3);
    if T < D
        Pop = false(N, D);
        for i = 1 : N
            k = randperm(T, 1);
            j = randperm(D, k);
            Pop(i, j) = 1;
        end
        Pop = Problem.Evaluation(Pop);
    else
        Pop = Problem.Initialization();
    end
end