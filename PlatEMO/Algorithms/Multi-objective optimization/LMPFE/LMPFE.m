classdef LMPFE < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% Evolutionary algorithm with local model based Pareto front estimation
% fPFE  --- 0.1 --- Frequency of employing generic front modeling
% K     --- 5  --- Number of subregions

%------------------------------- Reference --------------------------------
% Y. Tian, L. Si, X. Zhang, K. C. Tan, and Y. Jin, Local model based Pareto
% front estimation for multi-objective optimization, IEEE Transactions on
% Systems, Man, and Cybernetics: Systems, 2023, 53(1): 623-634.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm, Problem)
            %% Parameter setting
            [fPFE,K] = Algorithm.ParameterSet(0.1,5);

            %% Generate random population
            Population = Problem.Initialization();
            FrontNo    = NDSort(Population.objs,inf);

            %% Generate K subregions
            [Center,R] = adaptiveDivision(Population.objs,K);
            %% Calculate the intersection point on each subregion
            % Initialze parameter P
            P = ones(K,Problem.M);

            % Calculate the fitness of each solution
            [App,Dis] = subFitness(Population.objs,P,Center,R);
            Dis       = sort(Dis,2); 
            Crowd     = Dis(:,1) + 0.1*Dis(:,2);

            theta     = 0.8;
            preApp    = mean(App);
            preCrowd  = mean(Crowd);
            %% Optimization
            while Algorithm.NotTerminated(Population)
                % Mating
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-theta*Crowd-(1-theta)*App);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                % Generic front modeling
                if ~mod(ceil(Problem.FE/Problem.N),ceil(fPFE*ceil(Problem.maxFE/Problem.N))) || fPFE == 0
                    % Update subregions
                    [Center,R] = adaptiveDivision(Population.objs,K);
                    % PF modeling
                    P = subGFM(Population.objs,Center,R,FrontNo);
                end
                [Population,FrontNo,App,Crowd] = EnvironmentalSelection([Population,Offspring],P,theta,Problem.N,Center,R);
                % Update theta
                [theta,preApp,preCrowd] = UpdateTheta(preApp,preCrowd,App,Crowd);
            end
        end
    end
end