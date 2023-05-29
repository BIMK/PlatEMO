classdef PCSAEA < ALGORITHM
% <multi/many> <real> <expensive>
% Pairwise comparison based surrogate-assisted evolutionary algorithm
% delta ---  0.8 --- Threshold of reliability measurement
% gmax  --- 3000 --- Number of solutions evaluated by surrogate model

%------------------------------- Reference --------------------------------
% Y. Tian, J. Hu, C. He, H. Ma, L. Zhang, and X. Zhang, A pairwise
% comparison based surrogate-assisted evolutionary algorithm for expensive
% multi-objective optimization, Swarm and Evolutionary Computation, 2023.
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
            %% Parameter setting
            [delta,gmax] = Algorithm.ParameterSet(0.8,3000);

            %% Initalize the population by Latin hypercube sampling
            N          = max(11*Problem.D-1,Problem.N);
            PopDec     = UniformPoint(N,Problem.D,'Latin');
            Population = Problem.Evaluation(repmat(Problem.upper-Problem.lower,N,1).*PopDec+repmat(Problem.lower,N,1));
            Arc        = Population;
            t = 1;
        
            %% Optimization
            while Algorithm.NotTerminated(Arc)
                % Select a balance sample set by a new fitness
                [Input,Output,Pa,Pmid] = CalFitnessPC(Population.objs,Population.decs,(Problem.FE/Problem.maxFE));   
                % Data process
                [TrainIn,~,TestIn,TestOut] = DataProcess(Input,Output);
                % Construct and update the FNN£¬global classify surrogate model
                net = RBFNNPC(0.1925);             
                net.train(TrainIn,Problem.D);              

                % Error rates calculation
                TestPre   = net.lastpredict(TestIn,Problem.D,Pmid,1);
                % New and suitble reliability selection
                validIndex = TestPre~=1.5;
                Error1 = sum(TestOut(validIndex)==TestPre(validIndex))/length(TestOut);   
                Error2 = sum(TestOut(validIndex)~=TestPre(validIndex))/length(TestOut);              

                % Surrogate-assisted selection and update the population
                Next = SurrogateAssistedSelectionPC(Problem,net,Error1,Error2,Population.decs,gmax,Pa,Problem.D,0,delta);
                if ~isempty(Next)
                    Arc = [Arc,Problem.Evaluation(Next)];
                end

                Population = EnvironmentalSelection(Arc,Problem.N);                
                t = t+1;
            end
        end
    end
end