classdef MOEADDAE < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained>
% MOEA/D with detect-and-escape strategy

%------------------------------- Reference --------------------------------
% Q. Zhu, Q. Zhang, and Q. Lin, A constrained multi-objective evolutionary
% algorithm with detect-and-escape strategy, IEEE Transactions on
% Evolutionary Computation, 2020, 24(5): 938-947.
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
            %% Generate the weight vectors
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            % Size of neighborhood
            T = ceil(Problem.N/5);

            %% Detect the neighbours of each solution
            B = pdist2(W,W);
            [~,B] = sort(B,2);
            B = B(:,1:T);

            %% Generate random population
            Population  = Problem.Initialization();
            z           = min(Population.objs,[],1);
            Pi          = ones(Problem.N,1);
            CP_old      = sum(sum(max(Population.cons,0),2));
            avg_fit     = sum(max(abs((Population.objs-repmat(z,Problem.N,1)).*W),[],2))/Problem.N;
            epsilon_max = max(sum(max(Population.cons,0),2));
            epsilon     = epsilon_max*size(Population.cons,2);
            sigma_min   = 1-(1/(epsilon+1.0e-10))^(3/ceil(Problem.maxFE/Problem.N));
            current     = 0;
            A     = ArchiveUpdate(Population,Problem.N);
            gen   = 0;    
            CV    = sum(max(Population.cons,0),2);
            fr    = length(find(CV<=0))/Problem.N;
            sigma = max(sigma_min,fr);

            %% Optimization
            while Algorithm.NotTerminated(A)
                Q = [];
                for subgeneration = 1 : 5
                    Bounday = find(sum(W<1e-3,2)==Problem.M-1)';
                    I = [Bounday,TournamentSelection(10,floor(Problem.N/5)-length(Bounday),-Pi)];
                    for j = 1 : length(I)
                        i = I(j);
                        % Choose the parents
                        if rand < 0.9
                            P = B(i,randperm(size(B,2)));
                        else
                            P = randperm(Problem.N);
                        end
                        Offspring = OperatorGAhalf(Problem,Population(P(1:2)));
                        z = min([z;Offspring.objs],[],1);
                        [Population,Pi] = UpdatePop(Population,P,Offspring,epsilon,z,W,sigma,Pi,avg_fit);
                        if current == 1
                            tA = UpdatetA(tA,P,Offspring,z,W,avg_fit);
                        end   
                        if sum(max(Offspring.cons,0),2) == 0
                            Q = [Q Offspring];
                        end
                    end 
                end
                gen   = gen+1;
                CV    = sum(max(Population.cons,0),2);
                fr    = length(find(CV<=0))/Problem.N;
                sigma = max(sigma_min,fr);
                if current==0 || current==2
                    if fr < 0.95
                        epsilon = (1-sigma)*epsilon;
                    else
                        if current == 0
                           current = 1;
                           epsilon = 1e+30;
                           tA      = Population;
                           gen     = 1;
                        elseif current == 2
                           epsilon   = epsilon_max;
                           sigma_min = 1-(1/(epsilon+1.0e-10))^(3/ceil(Problem.maxFE/Problem.N));
                        end 
                    end 
                end
                if mod(gen,10) == 0
                    CP      = sum(sum(max(Population.cons,0)));
                    ROC     = abs(CP-CP_old)/(CP+1e-10);
                    CP_old  = CP;
                    avg_fit = sum(max(abs((Population.objs-repmat(z,Problem.N,1)).*W),[],2))/Problem.N;
                    if current == 0
                        if (ROC>0&&ROC<1e-5) && (CP>0.1*epsilon_max)
                            current = 1;
                            epsilon = 1e+30;
                            gen     = 1;
                            tA      = Population;
                        end 
                    elseif current == 1
                        if ROC>0 && ROC<1e-5
                           current     = 2;
                           Population  = tA;
                           epsilon_max = max(sum(max(Population.cons,0),2));
                           epsilon     = epsilon_max;
                           sigma_min   = 1-(1/(epsilon+1.0e-10))^(3/ceil(Problem.maxFE/Problem.N));
                           z           = min(Population.objs,[],1);
                        end
                    end
                end
                if size(Q,2) > 0
                   A = ArchiveUpdate([A Q],Problem.N);
                end
            end
        end
    end
end