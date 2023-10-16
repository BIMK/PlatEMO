classdef ESBCEO < ALGORITHM
% <multi> <real> <expensive>
% Bayesian co-evolutionary optimization based entropy search
% Ke    ---   5 --- The number of function evaluations at each generation
% delta --- 0.9 --- The probability of choosing parents locally
% nr    ---   2 --- Maximum number of solutions replaced by each offspring
% L1    ---  80 --- The maximal number of points used for building a local model
% L2    ---  20 --- The maximal number of points used for building a local model

%------------------------------- Reference --------------------------------
% H. Bian, J. Tian, J. Yu, and H. Yu, Bayesian co-evolutionary optimization
% based entropy search for high-dimensional many-objective optimization,
% Knowledge-Based Systems, 2023, 274: 110630.
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
            [Ke,delta,nr,L1,L2] = Algorithm.ParameterSet(5,0.9,2,80,20);

            %% Generate random population
            NI = 100+Problem.D/10;
            P  = UniformPoint(NI,Problem.D,'Latin');
            Population = Problem.Evaluation(repmat(Problem.upper-Problem.lower,NI,1).*P+repmat(Problem.lower,NI,1));
            L1 = min(L1,length(Population));

            %% Optimization
            while Algorithm.NotTerminated(Population)
                PopDec     = EGOSelect(Problem,Population,L1,L2,Ke,delta,nr);
                Offspring  = Problem.Evaluation(PopDec);
                Population = [Population,Offspring];
            end
        end
    end
end