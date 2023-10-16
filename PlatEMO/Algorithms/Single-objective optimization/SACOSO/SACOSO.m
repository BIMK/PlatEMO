classdef SACOSO < ALGORITHM
% <single> <real/integer> <large/none> <expensive>
% Surrogate-assisted cooperative swarm optimization
% NFES ---  30 --- Population size of FES-assisted PSO
% NRBF --- 200 --- Population size of RBFNN-assisted PSO

%------------------------------- Reference --------------------------------
% C. Sun, Y. Jin, R. Cheng, J. Ding, and J. Zeng, Surrogate-assisted
% cooperative swarm optimization of high-dimensional expensive problems,
% IEEE Transactions on Evolutionary Computation, 2017, 21(4): 644-660.
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
            assert(~isempty(ver('nnet')),'The execution of SA-COSO requires the Deep Learning Toolbox.');
            
            %% Parameter setting
            [NFES,NRBF] = Algorithm.ParameterSet(30,200);
            
            %% Generate the random population
            N          = NFES + NRBF;
            PopDec     = UniformPoint(N,Problem.D,'Latin');
            Population = Problem.Evaluation(repmat(Problem.upper-Problem.lower,N,1).*PopDec+repmat(Problem.lower,N,1));
            MaxNode    = 8;
            NDB        = MaxNode*Problem.D+10;
            BU         = Problem.upper;
            BD         = Problem.lower;
            
            %% Initialize each swarm
            % FES-assisted swarm
            PosFES   = Population(1:NFES).decs;
            ObjFES   = Population(1:NFES).objs;
            SwarmFES = [PosFES,ObjFES];
            PbestFES = SwarmFES;
            PBEval   = true(N,1);
            notEval  = false(NFES,1);
            notEST   = false(NFES,1);
            [~,best] = min(SwarmFES(:,end));
            GbestFES = SwarmFES(best,:);
            VelFES   = repmat(BU-BD,NFES,1).*rand(NFES,1)+repmat(BD,NFES,1);
            % RBF-assisted SL-PSO
            PosRBF   = Population(NFES+1:NFES+NRBF).decs;
            ObjRBF   = Population(NFES+1:NFES+NRBF).objs;
            SwarmRBF = [PosRBF,ObjRBF];
            [~,best] = min(SwarmRBF(:,end));
            GbestRBF = SwarmRBF(best,:);
            DeltaRBF = repmat(BU-BD,NRBF,1).*rand(NRBF,1)+repmat(BD,NRBF,1);
            
            %% Optimization
            SwarmFESt = SwarmFES;
            SpreadSum = 0;
            Archive   = Population;
            iter      = 1;
            while Algorithm.NotTerminated(Population)
                % Determine global best solution
                if GbestFES(:,end) < GbestRBF(:,end)
                    Gbest = GbestFES;
                else
                    Gbest = GbestRBF;
                end

                % Paramaters of RBFNN
                SpreadSum = SpreadSum + sqrt(sum((max(Population.decs,[],1)-min(Population.decs,[],1)).^2));
                Spread    = SpreadSum/iter;
                 
                % Build RBFNN surrogate model
                net = srgtsnewrb(Population.decs',Population.objs',0.1,Spread,MaxNode,1,'off');
                
                % FES-assisted PSO
                SwarmFESt1 = SwarmFESt;
                SwarmFESt  = SwarmFES;
                [tArchive,SwarmFES,VelFES,PbestFES,GbestFES,notEval,notEST,PBEval] = FESPSO(iter,net,SwarmFES,VelFES,PbestFES,GbestFES,GbestRBF,Problem,notEval,notEST,PBEval,SwarmFESt1,SwarmFESt);
                
                % Select demonstrators
                Select = randperm(length(Population),NRBF);
                Demons = [SwarmRBF;Population(Select).decs,Population(Select).objs];
                % RBF-assisted SL-PSO
                [SwarmRBF,DeltaRBF] = RBFOperator(net,Demons,SwarmRBF,DeltaRBF,Gbest,Problem);
                [tArchive,GbestRBF] = UpdateRBF(Problem,tArchive,SwarmRBF,GbestRBF);
                
                % Updating archive DB
                Population = [Population,tArchive];
                Archive    = UpdateArchive(Archive,tArchive,SwarmRBF,NDB);
                iter       = iter + 1;
            end           
        end
    end
end