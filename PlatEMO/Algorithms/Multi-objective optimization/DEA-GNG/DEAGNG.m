classdef DEAGNG < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% Decomposition based evolutionary algorithm guided by growing neural gas
% aph ---   0.1 --- Parameter alpha
% eps --- 0.314 --- Parameter epsilon

%------------------------------- Reference --------------------------------
% Y. Liu, H. Ishibuchi, N. Masuyama, and Y. Nojima, Adapting reference
% vectors and scalarizing functions by growing neural gas to handle
% irregular Pareto fronts. IEEE Transactions on Evolutionary Computation,
% 2020, 24(3): 439-453.
%--------------------------------------------------------------------------
% Copyright Yiping Liu
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [aph,eps] = Algorithm.ParameterSet(0.1,0.314); 
                           
            %% DEA Initialization
            [Ru,Problem.N] = UniformPoint(Problem.N,Problem.M);	% Uniform reference vectors
            Population     = Problem.Initialization();          % Random population
            Zmin           = min(Population.objs,[],1);         % Ideal Point  
            AS             = [];                                % Input Signal Archive
            Ruq            = Ru;                                % Refernce vectors in Ru for selection
            [FrontNo,~]    = NDSort(Population.objs,Problem.N); % Fitness for the first mating selection
            crd            = zeros(1,Problem.N);                % Fitness for the first mating selection   
            MaxGen         = ceil(Problem.maxFE/Problem.N);    	% Maximum Generation
            
            %% GNG Initialization
            ArchiveSize = Problem.M*Problem.N;	% Size of Input Signal Archive
            NoG = aph*MaxGen;                   % Number of generations of Not Training GNG
            GNGnet.maxIter = 1;                 % Number of iterations to train GNG per Generation  
            GNGnet.maxAge = Problem.N;          % Maximum cluster age     
            GNGnet.maxNode = Problem.N;         % Max number of nodes
            GNGnet.lambda = 0.2*Problem.N;      % Cycle for topology reconstruction   
            GNGnet.hp = [];                     % Hit point of node 
            GNGnet.maxHP = 2*ArchiveSize;       % Max HP of node 
            GNGnet.Node = [];                   % Node
            GNGnet.NodeS = [];                  % Expanded node 
            GNGnet.NodeP = [];                  % Node mapped to hyperplane 
            GNGnet.Err = [];                    % Error 
            GNGnet.edge = zeros(2,2);           % Edge between nodes 
            GNGnet.age = zeros(2,2);            % Age of edge 
            GNGnet.epsilon_a = 0.2;             % Learning coefficient
            GNGnet.epsilon_nb = 0.01;           % Learning coefficient of neighbor
            GNGnet.alpha = 0.5;                 % Nodes r1max and r2max error reduction constant
            GNGnet.delta = 0.9;                 % Error reduction coefficient 

            %% Optimization
            while Algorithm.NotTerminated(Population)        
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,crd);
                Offspring  = OperatorGA(Problem,Population(MatingPool));       
                Zmin       = min([Zmin;Offspring.objs],[],1);

                %% GNG-based adaptation
                if ceil(Problem.FE/Problem.N) <= MaxGen - NoG   
                    % Input Signal Archive Update
                    AS = ArchiveUpdate([AS;Offspring.objs],ArchiveSize,Ruq,GNGnet.NodeS,Zmin);
                    nAS = length(AS);      
                    % GNG Update (and Algorithm 3)
                    GNGnet.maxNode = min(Problem.N,floor(nAS/2)); % paramter reset 
                    GNGnet.maxHP = 2*nAS; % paramter reset
                    GNGnet = GNGUpdate(AS,GNGnet);
                    % Reference Vector Adaptation (Algorithm 4) 
                    if size(GNGnet.NodeS,1)>2
                        [Ruq,GNGnet] = ReferenceCombination(Ru,GNGnet);
                    end
                    % Scalarizing Function Adaptation
                    theta = TunePBI(GNGnet,eps); % Tune theta in PBI function                    
                end

               %% Environmental Selection                 
               [Population,FrontNo,crd] = ESelection([Population,Offspring],Problem.N,Ruq,GNGnet.NodeS,theta,Zmin);
            end
        end
    end
end