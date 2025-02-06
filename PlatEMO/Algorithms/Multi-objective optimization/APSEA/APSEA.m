classdef APSEA < ALGORITHM
% <2025> <multi> <real/integer/label/binary/permutation> <constrained>
% Adaptive population sizing based evolutionary algorithm
% alpha --- 0.05 --- Minimum ratio of feasible solution in Population1
% beta  --- 0.05 --- Minimum change rate of objectives in Population1
% cp    ---    5 --- Decrease trend of the dynamic constraint boundary

%------------------------------- Reference --------------------------------
% Y. Tian, R. Wang, Y. Zhang, and X. Zhang. Adaptive population sizing for
% multi-population based constrained multi-objective optimization.
% Neurocomputing, 2025: 129296.
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
            [alpha,beta,cp] = Algorithm.ParameterSet(0.05,0.05,5);
            Population1  = Problem.Initialization();
            Population2  = Problem.Initialization();               
            Fitness1     = CalFitness(Population1.objs,Population1.cons);
            Fitness2     = CalFitness(Population2.objs);  
            last_gen     = 20;
            ideal_points = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);
            nadir_points = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);
            epsilon0     = max(sum(max(0,Population2.cons),2));
            if epsilon0 == 0
                epsilon0 = 1;
            end
            gen = 1;

            %% Optimization
            while Algorithm.NotTerminated(Population1)
                FR = sum(sum(max(Population1.cons,0),2)==0)/length(Population1); % the ratio of feasible soluitons in Population1
                ideal_points(gen,:) = min(Population1.objs,[],1);
                nadir_points(gen,:) = max(Population1.objs,[],1);
                if gen > last_gen
                    max_change = calc_maxchange(ideal_points,nadir_points,gen,last_gen);
                end
                if FR <= alpha || gen <= last_gen
                    %% update the size of Population2
                    N = max(ceil(Problem.N/2*(1-log2(1+FR))+Problem.N/2*(1-log2(1+Problem.FE/Problem.maxFE))),2); 
                    %% reproduction and selection
                    MatingPool1 = TournamentSelection(2,Problem.N,Fitness1);
                    MatingPool2 = TournamentSelection(2,N,Fitness2);
                    Offspring1  = OperatorGA(Problem,Population1(MatingPool1));
                    Offspring2  = OperatorGA(Problem,Population2(MatingPool2));
                    [Population1, Fitness1] = EnvironmentalSelection_CDP([Population1,Offspring1,Offspring2],Problem.N);
                    [Population2, Fitness2] = EnvironmentalSelection_noConstrained([Population2,Offspring1,Offspring2,],N);
                elseif max_change > beta
                    %% update the epsilon-constraint boundary
                    epsilon = ReduceBoundary(epsilon0,gen,ceil(Problem.maxFE/Problem.N)-1,cp);
                    %% update the size of Population2
                    N = max(ceil(Problem.N/2*(1-log2(1+FR))+Problem.N/2*(1-log2(1+Problem.FE/Problem.maxFE))),2); 
                    %% reproduction and selection
                    MatingPool1 = TournamentSelection(2,Problem.N,Fitness1);
                    MatingPool2 = TournamentSelection(2,N,Fitness2);
                    Offspring1  = OperatorGA(Problem,Population1(MatingPool1));
                    Offspring2  = OperatorGA(Problem,Population2(MatingPool2));
                    [Population1, Fitness1] = EnvironmentalSelection_CDP([Population1,Offspring1,Offspring2],Problem.N);
                    [Population2, Fitness2] = EnvironmentalSelection_Epsilon([Population2,Offspring1,Offspring2],N,epsilon);  
                else
                    MatingPool1 = TournamentSelection(2,Problem.N,Fitness1);
                    Offspring1  = OperatorGA(Problem,Population1(MatingPool1));
                    [Population1, Fitness1] = EnvironmentalSelection_CDP([Population1,Offspring1],Problem.N);
                end
                gen = gen+1;
            end
        end
    end                
end

function max_change = calc_maxchange(ideal_points,nadir_points,gen,last_gen)
% Calculate the maximum rate of change of objectives

    delta_value = 1e-6 * ones(1,size(ideal_points,2));
    rz  = abs((ideal_points(gen,:) - ideal_points(gen - last_gen + 1,:)) ./ max(ideal_points(gen - last_gen + 1,:),delta_value));
    nrz = abs((nadir_points(gen,:) - nadir_points(gen - last_gen + 1,:)) ./ max(nadir_points(gen - last_gen + 1,:),delta_value));
    max_change = max([rz, nrz]);
end

function epsn = ReduceBoundary(eF, k, MaxK, cp)
% Update the epsilon constraint boundary

    z        = 1e-8;
    Nearzero = 1e-15;
    B        = MaxK./power(log((eF + z)./z), 1.0./cp);
    B(B==0)  = B(B==0) + Nearzero;
    f        = eF.* exp( -(k./B).^cp );
    tmp      = find(abs(f-z) < Nearzero);
    f(tmp)   = f(tmp).*0 + z;
    epsn     = f - z;
    epsn(epsn<=0) = 0;
end
