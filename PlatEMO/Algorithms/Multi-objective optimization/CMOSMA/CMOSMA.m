classdef CMOSMA < ALGORITHM
% <multi/many> <real/integer> <constrained>
% Constrained multi-objective evolutionary algorithm with self-organizing map

%------------------------------- Reference --------------------------------
% C. He, M. Li, C. Zhang, H. Chen, P. Zhong, Z. Li, and J. Li, A
% self-organizing map approach for constrained multi-objective optimization
% problems, Complex & Intelligent Systems, 2022, 8: 5355-5375.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Chao He

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [D,tau0,H] = Algorithm.ParameterSet(repmat(ceil(Problem.N.^(1/(Problem.M-1))),1,Problem.M-1),0.7,5);
            D1         = D;
            Problem.N  = prod(D);
            sigma0     = sqrt(sum(D.^2)/(Problem.M-1))/2;

            %% Generate random population
            FP = Problem.Initialization();
            AP = Problem.Initialization();

            %% Initialize the SOM
            % Training set
            S = FP.decs;
            % Weight vector of each neuron
            W = S;
            [LDis,B] = Initialize_SOM(S,D,H);
            % Training set
            S2 = AP.decs;
            % Weight vector of each neuron
            W2 = S2;
            [LDis2,B2] = Initialize_SOM(S2,D1,H);
            
            %% Optimization
            while Algorithm.NotTerminated(FP)
                % Update SOM1
                W = UpdateSOM(S,W,Problem.FE,Problem.maxFE,LDis,sigma0,tau0); 
                % Update SOM2
                W2 = UpdateSOM(S2,W2,Problem.FE,Problem.maxFE,LDis2,sigma0,tau0);

                % Associate each solution with a neuron
                XU  = Associate(FP,W,Problem.N);
                XU2 = Associate(AP,W2,Problem.N); 

                %Construct  matingPool  
                [MatingPool1] = MatingPool(XU,Problem.N,B);
                [MatingPool2] = MatingPool(XU2,Problem.N,B2);

                % Evolution
                A1 = FP.decs;
                Offspring1 = OperatorGA(Problem,[FP(XU),FP(MatingPool1)]);%GA
                A2         = AP.decs;
                Offspring2 = OperatorGA(Problem,[AP(XU2),AP(MatingPool2)]);%GA
                %EnvironmentalSelection
                FP = EnvironmentalSelection([FP,Offspring1,Offspring2],Problem.N,true);
                AP = EnvironmentalSelection([AP,Offspring1,Offspring2],Problem.N,false);
                % Update the training set
                S  = setdiff(FP.decs,A1,'rows');
                S2 = setdiff(AP.decs,A2,'rows');
            end
        end
    end
end