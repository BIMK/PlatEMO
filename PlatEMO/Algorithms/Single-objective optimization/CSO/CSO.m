classdef CSO < ALGORITHM
% <single> <real/integer> <large/none> <constrained/none>
% Competitive swarm optimizer
% phi --- 0.1 --- Social factor

%------------------------------- Reference --------------------------------
% R. Cheng and Y. Jin, A competitive swarm optimizer for large scale
% optimization, IEEE Transactions on Cybernetics, 2014, 45(2): 191-204.
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
            phi = Algorithm.ParameterSet(0.1);
            
            %% Generate random population
            Population = Problem.Initialization();
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                % Determine the losers and winners
                rank    = randperm(Problem.N);
                loser   = rank(1:end/2);
                winner  = rank(end/2+1:end);
                replace = FitnessSingle(Population(loser)) < FitnessSingle(Population(winner));
                temp            = loser(replace);
                loser(replace)  = winner(replace);
                winner(replace) = temp;
                % Update the losers by learning from the winners
                LoserDec  = Population(loser).decs;
                WinnerDec = Population(winner).decs;
                LoserVel  = Population(loser).adds(zeros(size(LoserDec)));
                R1 = rand(Problem.N/2,Problem.D);
                R2 = rand(Problem.N/2,Problem.D);
                R3 = rand(Problem.N/2,Problem.D);
                LoserVel = R1.*LoserVel + R2.*(WinnerDec-LoserDec) + phi.*R3.*(repmat(mean(Population.decs,1),Problem.N/2,1)-LoserDec);
                LoserDec = LoserDec + LoserVel;
                Population(loser) = Problem.Evaluation(LoserDec,LoserVel);
            end
        end
    end
end