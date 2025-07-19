classdef DMMOEA < ALGORITHM
% <2025> <multi> <real/integer/binary> <large/none> <constrained/none> <dynamic> <sparse>
% Dual model based multi-objective evolutionary algorithm

%------------------------------- Reference --------------------------------
% P. Zhang, R. Zhang, Y. Tian, K. C. Tan, and X. Zhang. A dual model-based
% evolutionary framework for dynamic large-scale sparse multiobjective
% optimization. Swarm and Evolutionary Computation, 2025.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
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
            ChangeCount = 0;
            AllPop      = [];
            DecSource   = cell(1,200);
            MaskSource  = cell(1,200);
            NDS         = cell(1,200);
            % Reset the number of saved populations (only for dynamic optimization)
            Algorithm.save = sign(Algorithm.save)*inf;
            
            %% Generate random population
            TDec    = [];
            TMask   = [];
            TempPop = [];
            Fitness = zeros(1,Problem.D);
            for i = 1 : 1+4*any(Problem.encoding~=4)
                Dec = unifrnd(repmat(Problem.lower,Problem.D,1),repmat(Problem.upper,Problem.D,1));
                Dec(:,Problem.encoding==4) = 1;
                Mask       = eye(Problem.D);
                Population = Problem.Evaluation(Dec.*Mask);
                TDec       = [TDec;Dec];
                TMask      = [TMask;Mask];
                TempPop    = [TempPop,Population];
                Fitness    = Fitness + NDSort([Population.objs,Population.cons],inf);
            end
            % Generate initial population
            Dec = unifrnd(repmat(Problem.lower,Problem.N,1),repmat(Problem.upper,Problem.N,1));
            Dec(:,Problem.encoding==4) = 1;
            Mask = false(Problem.N,Problem.D);
            for i = 1 : Problem.N
                Mask(i,TournamentSelection(2,ceil(rand*Problem.D),Fitness)) = 1;
            end
            Population    = Problem.Evaluation(Dec.*Mask);        
            DecSource{1}  = Dec;
            MaskSource{1} = Mask;
            NDS{1}        = Population;
            [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection([Population,TempPop],[Dec;TDec],[Mask;TMask],Problem.N);   

            %% Optimization
            while Algorithm.NotTerminated(Population)
                if Changed(Problem,Population)
                    PopDec      = Dec(FrontNo==1,:);
                    PopMask     = Mask(FrontNo==1,:);
                    Pop         = Problem.Evaluation(PopDec.*PopMask);
                    ChangeCount = ChangeCount+1;
                    DecSource{ChangeCount+1}  = Dec;
                    MaskSource{ChangeCount+1} = Mask;
                    NDS{ChangeCount+1} = Population;
                    AllPop = [AllPop,Population];
                    % React to the change
                    [Population,Dec,Mask] = Prediction(Problem,ChangeCount,DecSource,MaskSource);
                    [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Pop],[Dec;PopDec],[Mask;PopMask],Problem.N);                           
                end 
                
                MatingPool       = TournamentSelection(2,2*Problem.N,FrontNo,-CrowdDis);
                [OffDec,OffMask] = Operator(Problem,Dec(MatingPool,:),Mask(MatingPool,:),Fitness);
                Offspring        = Problem.Evaluation(OffDec.*OffMask);  
                [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Problem.N);

                if Problem.FE >= Problem.maxFE
                    % Return all populations
                    Population = [AllPop,Population];
                    [~,rank]   = sort(Population.adds(zeros(length(Population),1)));
                    Population = Population(rank);
                end 
            end
        end
    end
end