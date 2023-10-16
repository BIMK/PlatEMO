classdef MOEAPSL < ALGORITHM
% <multi> <real/integer/binary> <large/none> <constrained/none> <sparse>
% Multi-objective evolutionary algorithm based on Pareto optimal subspace
% learning

%------------------------------- Reference --------------------------------
% Y. Tian, C. Lu, X. Zhang, K. C. Tan, and Y. Jin, Solving large-scale
% multi-objective optimization problems with sparse optimal solutions via
% unsupervised neural networks, IEEE Transactions on Cybernetics, 2021,
% 51(6): 3115-3128.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Population initialization
            P   = UniformPoint(Problem.N,Problem.D,'Latin');
            Dec = P.*repmat(Problem.upper-Problem.lower,Problem.N,1) + repmat(Problem.lower,Problem.N,1);
            Dec(:,Problem.encoding==4) = 1;
            Mask = UniformPoint(Problem.N,Problem.D,'Latin') > 0.5;
            Population = Problem.Evaluation(Dec.*Mask);
            [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Dec,Mask,Problem.N,0,0);

            %% Optimization
            rho = 0.5;
            while Algorithm.NotTerminated(Population)
                Site = rho > rand(1,ceil(Problem.N/2));
                if any(Site)
                    [rbm,dae,allZero,allOne] = ModelTraining(Mask,Dec,any(Problem.encoding~=4));
                else
                    [rbm,dae,allZero,allOne] = deal([]);
                end
                MatingPool = TournamentSelection(2,ceil(Problem.N/2)*2,FrontNo,-CrowdDis);
                [OffDec,OffMask] = Operator(Problem,Dec(MatingPool,:),Mask(MatingPool,:),rbm,dae,Site,allZero,allOne);
                Offspring = Problem.Evaluation(OffDec.*OffMask);
                [Population,Dec,Mask,FrontNo,CrowdDis,sRatio] = EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Problem.N,length(Population),2*sum(Site));
                rho = (rho+sRatio)/2;
            end
        end
    end
end