classdef EGO < ALGORITHM
% <single> <real/integer> <expensive>
% Efficient global optimization
% IFEs --- 10000 --- Internal GA evals per iteration

%------------------------------- Reference --------------------------------
% D. R. Jones, M. Schonlau, and W. J. Welch, Efficient global optimization
% of expensive black-box functions, Journal of Global optimization, 1998,
% 13(4): 455-492.
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
            IFEs = Algorithm.ParameterSet(10000);
            
            %% Generate the random population
            N          = 10*Problem.D;
            PopDec     = UniformPoint(N,Problem.D,'Latin');
            Population = Problem.Evaluation(repmat(Problem.upper-Problem.lower,N,1).*PopDec+repmat(Problem.lower,N,1));
            theta      = 10.*ones(1,Problem.D);
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                % Select initial population
                [N,D] = size(Population.decs);
                if N > 11*D-1
                    [~,index] = sort(Population.objs);
                    Next      = index(1:11*D-1);
                else
                    Next = true(N,1);
                end
                PDec = Population(Next).decs;
                PObj = Population(Next).objs;
                
                % Surrogate-assisted prediction
                dmodel     = dacefit(Population.decs,Population.objs,'regpoly1','corrgauss',theta,1e-5.*ones(1,D),20.*ones(1,D));
                theta      = dmodel.theta;
                PopDec     = EvolEI(Problem,PDec,PObj,dmodel,IFEs);
                
                if checkExist(Population.decs,PopDec)
                    % Evaluate new candidate
                    Population = [Population,Problem.Evaluation(PopDec)];
                end
            end
        end
    end
end