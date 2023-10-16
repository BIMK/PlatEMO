classdef OSPNSDE < ALGORITHM
% <multi> <real/integer>
% Non-dominated sorting differential evolution with prediction in the objective space
% lambda --- 0.2 --- Hypervolume variation percent
% p      ---  50 --- Initial forecast horizon

%------------------------------- Reference --------------------------------
% E. Guerrero-Pena, A. F. R. Araujo, Multi-objective evolutionary
% algorithm with prediction in the objective space, Information Sciences,
% 2019, 501: 293-316.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Elaine Guerrero-Pena

    methods
        function main(Algorithm,Problem)
            assert(~isempty(ver('ident')),'The execution of OSP-NSDE requires the System Identification Toolbox.');
            
            %% Parameter setting
            [lambda,p] = Algorithm.ParameterSet(0.2,50);

            %% Generate random population
            Population = Problem.Initialization();
            FrontNo = NDSort(Population.objs,Population.cons,Problem.N);

            Nadyr_P = max(Population.objs);

            t_init = 2;  % Initial generation for regression
            t = 0;

            %% Optimization
            while Algorithm.NotTerminated(Population)
                t = t+1;
                Pt(:,:,t)     = Population.decs;
                Pt_apt(:,:,t) = Population.objs;
                Pt_Rank(:,t)  = FrontNo;
                Hyp(t,1) = HV(Population(FrontNo==1),Nadyr_P);

                if t==t_init
                    Nadyr_P = max(Pt_apt(:,:,t));
                end

                if t-t_init>=3 

                    alpha = Hyp(t_init) + lambda*Hyp(t_init);

                    % Criterion to trigger OSP procedure
                    if  t==10 || (Hyp(t) >= alpha && length(Pt_apt(FrontNo==1,:,t))<.9*Problem.N)

                        [Population, FrontNo] = OSP(t, Pt, Pt_apt, Pt_Rank,...
                            Problem, t_init, p);

                        % Update OSP criterion parameters and OSP parameter
                        t_init = t+1;
                        p = round(p/2);
                        if p==0
                            p=1;
                        end

                    else
                        MatingPool = [1:Problem.N,randi(Problem.N,1,2*Problem.N)];
                        Offspring = OperatorDE(Problem,Population(MatingPool(1:length(MatingPool)/3)),...
                            Population(MatingPool(length(MatingPool)/3+1:length(MatingPool)/3*2)),...
                            Population(MatingPool(length(MatingPool)/3*2+1:end)));              
                        [Population,FrontNo] = EnvironmentalSelection([Population,Offspring],Problem.N);
                    end
                else
                    MatingPool = [1:Problem.N,randi(Problem.N,1,2*Problem.N)];
                    Offspring = OperatorDE(Problem,Population(MatingPool(1:length(MatingPool)/3)),...
                        Population(MatingPool(length(MatingPool)/3+1:length(MatingPool)/3*2)),...
                        Population(MatingPool(length(MatingPool)/3*2+1:end)));               
                    [Population,FrontNo] = EnvironmentalSelection([Population,Offspring],Problem.N);
                end
            end
        end
    end
end