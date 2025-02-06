classdef LCMEA < ALGORITHM
% <2024> <multi> <real> <large/none> <constrained>
% Large-scale constrained multi-objective evolutionary algorithm

%------------------------------- Reference --------------------------------
% L. Si, X. Zhang, Y. Zhang, S. Yang, and Y. Tian. An efficient sampling
% approach to offspring generation for evolutionary large-scale constrained
% multi-objective optimization. IEEE Transactions on Emerging Topics in
% Computational Intelligence, 2024.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evoluationary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Set the default parameters
            % Initialize archive
            Archive = Problem.Initialization(2*Problem.N);

            % Environmental selection
            Env = EnvRL(Problem, Archive, floor((Problem.maxFE - Problem.FE) / Problem.N));

            % Parameter of relaxation method
            VAR0 = max(sum(max(0,Archive.cons),2));
            if VAR0 == 0
                VAR0 = 1;
            end
            X   = 0;
            cp  = (-log(VAR0)-6)/log(1-0.5);
            VAR = VAR0*(1-X)^cp;

            [Population, PopIndex] = Env.FormPop(Archive);
            PreAction = inf;
            while Algorithm.NotTerminated(Archive)
                % Generate offspring
                Offspring = ESP(Problem, Population);

                % Environmental selection         
                [Env, Population, Action] = Env.do([Population, Offspring], VAR, Archive);

                % Update Archive
                Archive(PopIndex) = Population;

                % Update Population
                if Action ~= PreAction
                    [Population, PopIndex] = Env.FormPop(Archive);
                end
                PreAction = Action;

                cp  = (-log(VAR0)-6)/log(1-0.5);
                VAR = VAR0*(1-X)^cp;
                if VAR < 1e-6
                    VAR = 0;
                end
                X = X+1/(Problem.maxFE/Problem.N);
            end
        end
    end
end