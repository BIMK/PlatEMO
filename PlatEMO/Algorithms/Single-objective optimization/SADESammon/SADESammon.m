classdef SADESammon < ALGORITHM
% <single> <real/integer> <expensive>
% Sammon mapping assisted differential evolution

%------------------------------- Reference --------------------------------
% G. Chen, K. Zhang, X. Xue, L. Zhang, J. Yao, H. Sun, L. Fan, and Y. Yang,
% Surrogate-assisted evolutionary algorithm with dimensionality reduction
% method for water flooding production optimization, Journal of Petroleum
% Science and Engineering, 2020, 185: 106633.
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
            [Mu,CR,LD,Omiga,Lammda] = Algorithm.ParameterSet(0.5,0.5,4,2,50);
            
            %% Generate the random population
            N          = 2*Problem.D;
            PopDec     = UniformPoint(N,Problem.D,'Latin');
            Population = Problem.Evaluation(repmat(Problem.upper-Problem.lower,N,1).*PopDec+repmat(Problem.lower,N,1));
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                % Generate Lammda new solution with DE
                Lammda = min(Lammda,length(Population));
                newDec = GenerateNew(Problem,Population,Lammda,Mu,CR);
                
                % Map data to low dimension using Sammon mapping
                trainDec = Population(end-N+1:end).decs;
                trainObj = Population(end-N+1:end).objs;
                sdrDec   = Sammon([trainDec;newDec],LD);
                
                % Build Kriging surrogate model
                try
                    dmodel = dacefit(sdrDec(1:N,:),trainObj,'regpoly1','corrgauss',0.002);
                    [yhat,predvar] = predictor(sdrDec(N+1:end,:), dmodel);
                catch
                    continue;
                end
                
                % Obtain LCB
                LCB = yhat - Omiga*predvar;
                
                % Evaluate
                [~,best] = min(LCB);
                Population = [Population,Problem.Evaluation(newDec(best,:))];
            end
        end
    end
end