classdef MOEADMRDL < ALGORITHM
% <multi> <real/integer>
% MOEA/D with maximum relative diversity loss
% gamma --- 20 --- Maximum allowable maximum relative diversity loss
% nmov  --- 10 --- Size of moving average

%------------------------------- Reference --------------------------------
% S. B. Gee, K. C. Tan, V. A. Shim, and N. R. Pal, Online diversity
% assessment in evolutionary multiobjective optimization: A geometrical
% perspective, IEEE Transactions on Evolutionary Computation, 2015, 19(4):
% 542-559.
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
            [gamma,nmov] = Algorithm.ParameterSet(20,10);

            %% Generate the weight vectors
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);	% Weight vectors
            T = ceil(Problem.N/10);                             % Neighbourhood size

            %% Detect the neighbours of each solution
            B = pdist2(W,W);
            [~,B] = sort(B,2);
            B = B(:,1:T)';                          % The neighbours of each weight vector

            %% Generate random population
            Population = Problem.Initialization();	% The initial population
            Z = min(Population.objs,[],1);          % The ideal point

            %% Optimization
            disC      = 20;     % The distribution index of SBX
            Pn        = 0;      % The standard deviation of Gaussian permutation
            allEgamma = [];     % The MRDL array storing all the past MRDL values
            while Algorithm.NotTerminated(Population)
                % Select parents
                MatingPool = [B(randi(T,1,Problem.N)+(0:Problem.N-1)*T),B(randi(T,1,Problem.N)+(0:Problem.N-1)*T)];
                % Generate offsprings
                OffDec    = OperatorGAhalf(Problem,Population(MatingPool).decs,{1,disC,1,20});
                OffDec    = OffDec + randn(size(OffDec))*Pn;
                Offspring = Problem.Evaluation(OffDec);
                % Update the ideal point
                Z          = min([Z;Offspring.objs],[],1);
                % Environmental selection
                [Population,Egamma] = EnvironmentalSelection(Population,Offspring,W,B',Z,gamma);
                % Operator adaption
                [disC,Pn,allEgamma] = OperatorAdaption(disC,Pn,allEgamma,Egamma,nmov);
            end
        end
    end
end