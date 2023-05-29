classdef SGECF < ALGORITHM
% <multi> <real/binary> <large/none> <constrained/none> <sparse>
% Sparsity-guided elitism co-evolutionary framework

%------------------------------- Reference --------------------------------
% C. Wu, Y. Tian, Y. Zhang, H. Jiang, and X. Zhang, A sparsity-guided
% elitism co-evolutionary framework for sparse large-scale multi-objective
% optimization, Proceedings of the IEEE Congress on Evolutionary
% Computation, 2023.
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
            % Calculate the fitness of each decision variable
            TDec    = [];
            TMask   = [];
            TempPop = [];
            Fitness = zeros(1,Problem.D);
            REAL    = ~strcmp(Problem.encoding,'binary');
            for i = 1 : 1+3*REAL
                if REAL
                    Dec = unifrnd(repmat(Problem.lower,Problem.D,1),repmat(Problem.upper,Problem.D,1));
                else
                    Dec = ones(Problem.D,Problem.D);
                end
                Mask       = eye(Problem.D);
                Population = Problem.Evaluation(Dec.*Mask);
                TDec       = [TDec;Dec];
                TMask      = [TMask;Mask];
                TempPop    = [TempPop,Population];
                Fitness    = Fitness + NDSort([Population.objs,Population.cons],inf);
            end
            % Calculate the theta
            if REAL
                VDec = unifrnd(repmat(Problem.lower,Problem.D,1),repmat(Problem.upper,Problem.D,1));
            else
                VDec = ones(Problem.D,Problem.D);
            end
            VMask = zeros(Problem.D,Problem.D);
            for i = 1 : Problem.D
                [~,rank1] = sort(Fitness);
                VMask(i,rank1(1:i)) = 1;
            end
            VPopulation = Problem.Evaluation(VDec.*VMask);
            [VPopulation,VDec,VMask,FrontNo,CrowdDis] = SPEA2_EnvironmentalSelection([VPopulation,TempPop],[VDec;TDec],[VMask;TMask],Problem.N);
            pop = VPopulation(find(FrontNo ==1));
            popMask     = VMask(find(FrontNo ==1),:);
            thetaAll    = sum(popMask,2)';
            thetaUnique = unique(thetaAll);
            if mod(size(thetaUnique,2),2)==0
                minOdd = max(1,size(thetaAll,2)-1);
            else
                minOdd = max(1,size(thetaAll,2)-2);
            end
            [~,theta] = kmeans(thetaAll',minOdd);
            thetamid  = median(theta);
            
            % Generate initial population
            Dec  = [];
            Mask = [];
            % Generate initial population
            if REAL
                Dec = unifrnd(repmat(Problem.lower,Problem.N,1),repmat(Problem.upper,Problem.N,1));
            else
                Dec = ones(Problem.N,Problem.D);
            end
            Mask = zeros(Problem.N,Problem.D);
            for i = 1 : Problem.N
                Mask(i,TournamentSelection(2,ceil(rand*Problem.D),Fitness)) = 1;
            end
            Population = Problem.Evaluation(Dec.*Mask);
            [Population,Dec,Mask,FrontNo,CrowdDis] = SPEA2_EnvironmentalSelection([VPopulation,Population],[VDec;Dec],[VMask;Mask],Problem.N);
            
            temptheta = sum(Mask,2)';
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                winIndex    = find(FrontNo ==1);
                winPopMask  = Mask(find(FrontNo ==1),:);
                winTheta    = sum(winPopMask,2)';
                thetaUnique = unique(winTheta);
                if mod(size(thetaUnique,2),2)==0
                    minOdd = max(1,size(winTheta,2)-1);
                else
                    minOdd = max(1,size(winTheta,2)-2);
                end
                [~,theta] = kmeans(winTheta',minOdd);
                thetamid  = median(theta);
                
                loseIndex   = find(FrontNo ~=1);
                losePopMask = Mask(find(FrontNo ~=1),:);
                loseTheta   = sum(losePopMask,2)';
                
                loseIndex1 = loseIndex(find(loseTheta<=thetamid));
                loseIndex2 = loseIndex(find(loseTheta>thetamid));
                
                [OffDec,OffMask] = Operator(Population,Dec,Mask,Fitness,winIndex,loseIndex1,loseIndex2,Problem,thetamid,REAL);
                Offspring        = Problem.Evaluation(OffDec.*OffMask);
                [Population,Dec,Mask,FrontNo,CrowdDis] = SPEA2_EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Problem.N);
                Theta = sum(Mask,2)';
            end
        end
    end
end