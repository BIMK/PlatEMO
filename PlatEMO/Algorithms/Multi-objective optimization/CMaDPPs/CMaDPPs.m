classdef CMaDPPs < ALGORITHM
% <2023> <multi/many> <real/integer/label/binary/permutation> <constrained/none>
% Constrained many-objective optimization with determinantal point processes

%------------------------------- Reference --------------------------------
% F. Ming, W. Gong, S. Li, L. Wang, and Z. Liao. Handling constrained
% many-objective optimization problems via determinantal point processes.
% Information Sciences, 2023, 643: 119260.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming (email: 20151000334@cug.edu.cn)

    methods
        function main(Algorithm,Problem)
            
            %% Parameter setting
            theta = Algorithm.ParameterSet(0);
            
            %% Generate random population
            Population = Problem.Initialization();
            CSA = Population;
            CV = sum(max(0,Population.cons),2);
            [z,znad] = deal(min(Population.objs),max(Population.objs));
            
            % parameters for epsilon
            CVmax = max(CV);
            G = ceil(Problem.maxFE/Problem.N);
            Tc = 0.8 * G;
            cp = 2;
            alpha = 0.95;
            tao = 0.05;
            epsilon = inf;
            
            % parameters for change
            last_gen         = 20;
            change_threshold = 1e-1;
            max_change       = 1;
            ideal_points     = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);
            nadir_points     = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);
            [Population,Archive] = EnvironmentalSelection(Population,[],[],Problem.N,z,znad,theta,CSA.objs,epsilon);       
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                gen        = ceil(Problem.FE/Problem.N);
                CV = sum(max(0,Population.cons),2);
                CV_max = max(CV);
                CVmax = max([CV_max,CVmax]);
                rf         = sum(CV <= 1e-6) / Problem.N;
                ideal_points(gen,:) = z;
                nadir_points(gen,:) = znad;
                
                % The maximumrate of change of ideal and nadir points rk is calculated.
                if gen >= last_gen
                    max_change = calc_maxchange(ideal_points,nadir_points,gen,last_gen);
                end
                
                % update epsilon before selection
                if max_change > change_threshold && gen < 0.4 * G
                    epsilon = inf;
                else
                    epsilon =  update_epsilon(tao,epsilon,CVmax,rf,alpha,gen,Tc,cp,Problem.M);
                end
                
                ParentC = MatingSelection([Population,Archive],CSA,Problem.N,z,znad);
                Offspring = OperatorGA(Problem,ParentC);
                
                z = min(z,min(Offspring.objs,[],1));
                CSA = UpdateCSA(CSA,Offspring,Problem.N,epsilon);
                
                [Population,Archive] = EnvironmentalSelection(Population,Offspring,Archive,Problem.N,z,znad,theta,CSA.objs,epsilon);
                znad = max(znad,max(Population.objs,[],1));
                if Problem.FE >= Problem.maxFE
                    Population = Archive;
                end
            end
        end
    end
end

% Calculate the Maximum Rate of Change
function max_change = calc_maxchange(ideal_points,nadir_points,gen,last_gen)
    delta_value = 1e-6 * ones(1,size(ideal_points,2));
    rz = abs((ideal_points(gen,:) - ideal_points(gen - last_gen + 1,:)) ./ max(ideal_points(gen - last_gen + 1,:),delta_value));
    nrz = abs((nadir_points(gen,:) - nadir_points(gen - last_gen + 1,:)) ./ max(nadir_points(gen - last_gen + 1,:),delta_value));
    max_change = max([rz, nrz]);
end

% Update the value of epsilon
function epsilon = update_epsilon(tao,epsilon_k,epsilon_0,rf,alpha,gen,Tc,cp,M)
    if gen > Tc
        epsilon = 0;
    else
        if rf < alpha
            epsilon = (1 - tao)^cp * epsilon_k;
        else
            epsilon = cp ^ M * epsilon_0 * ((1 - (gen / Tc)) ^ cp);
        end
    end
end