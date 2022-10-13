classdef TriMOEATAR < ALGORITHM
% <multi> <real/integer> <multimodal>
% Multi-modal MOEA using two-archive and recombination strategies
% p_con       --- 0.5  --- Probability of selecting parents from the convergence archive
% sigma_niche --- 0.1  --- Niche radius in the decision space
% eps_peak    --- 0.01 --- Accuracy level to detect peaks
% NR          --- 100  --- Number of refernece points
% NCA         --- 20   --- Number of sampling solutions in control variable analysis
% NIA         --- 6    --- Maximum number of tries required to judge the interaction

%------------------------------- Reference --------------------------------
% Y. Liu, G. G. Yen, and D. Gong, A multi-modal multi-objective
% evolutionary algorithm using two-archive and recombination strategies,
% IEEE Transactions on Evolutionary Computation, 2019, 23(4): 660-674.
%------------------------------- Copyright --------------------------------
% Copyright 2017-2018 Yiping Liu
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [p_con,sigma_niche,eps_peak,NR,NCA,NIA] = Algorithm.ParameterSet(0.5,0.1,0.01,100,20,6);  

            %% Decision Variable Analysis and Generate random population
            [Xic,Xre] = DecisionVariableAnalysis(Problem,NCA,NIA);

            %% Generate the reference points    
            R = UniformPoint(NR,Problem.M);

            %% Initialization
            Population = Problem.Initialization();
            NC = Problem.N;
            ND = Problem.N;
            AC = [];
            AD = [];
            Pc = floor(p_con.*Problem.N);
            Pd = Problem.N - Pc;

            %% Update Archives
            Z            = min(Population.objs,[],1);  
            [AC,RankC,~] = UpdateConvergenceArchive(AC,Population,NC,Z,Xic,sigma_niche,Problem);
            [AD,RankD]   = UpdateDiversityArchive(AD,Population,ND,R,Z,Xre,sigma_niche,Problem);

            %% Optimization
            while Algorithm.NotTerminated(AD)
                % Generate new Population based on two Archives
                MatingPoolC = TournamentSelection(2,Pc,RankC);
                MatingPoolD = TournamentSelection(2,Pd,RankD);
                Parents     = [AC(MatingPoolC),AD(MatingPoolD)];
                Parents     = Parents(randperm(Problem.N));
                Population  = OperatorGA(Problem,Parents);        
                % Update Archives
                Z             = min([Z;Population.objs],[],1);
                [AC,RankC,fS] = UpdateConvergenceArchive(AC,Population,NC,Z,Xic,sigma_niche,Problem);        
                [AD,RankD]    = UpdateDiversityArchive(AD,Population,ND,R,Z,Xre,sigma_niche,Problem);                
                % Recombination
                if Problem.FE >= Problem.maxFE && sum(Xic) > 0
                    FS = Recombination(AC,AD,Xic,Xre,eps_peak,fS,RankC);
                    AD = Problem.Evaluation(FS);
                end
            end
        end
    end
end