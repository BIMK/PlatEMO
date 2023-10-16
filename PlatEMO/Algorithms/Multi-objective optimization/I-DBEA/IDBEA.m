classdef IDBEA < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation> <constrained/none>
% Improved decomposition-based evolutionary algorithm

%------------------------------- Reference --------------------------------
% M. Asafuddoula, T. Ray, and R. Sarker, A decomposition-based evolutionary
% algorithm for many objective optimization, IEEE Transactions on
% Evolutionary Computation, 2015, 19(3): 445-460.
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
            %% Generate the reference points (vectors)
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            W = W./repmat(sqrt(sum(W.^2,2)),1,size(W,2));

            %% Generate random population
            Population = Problem.Initialization();
            z = min(Population.objs);
            a = Intercepts(Population.objs);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                % For each solution
                for i = 1 : Problem.N
                    % Generate an offspring
                    Offspring = OperatorGAhalf(Problem,Population([i,randi(Problem.N)]));

                    % Update the population
                    Feasible = all(Population.cons<=0,2);
                    if ~any(Feasible) || ~any(all(Population(Feasible).objs<=repmat(Offspring.obj,sum(Feasible),1),2))
                        % Calculate d1 and d2 values
                        List    = randperm(Problem.N);
                        nPopObj = (Population(List).objs-repmat(z,Problem.N,1))./repmat(a-z,Problem.N,1);
                        nOffObj = (Offspring.obj-z)./(a-z);
                        normP   = sqrt(sum(nPopObj.^2,2));
                        normO   = sqrt(sum(nOffObj.^2,2));
                        CosineP = sum(nPopObj.*W(List,:),2)./normP;
                        CosineO = sum(repmat(nOffObj,Problem.N,1).*W(List,:),2)./normO;
                        d1_old  = normP.*CosineP;
                        d1_new  = normO.*CosineO;
                        d2_old  = normP.*sqrt(1-CosineP.^2);
                        d2_new  = normO.*sqrt(1-CosineO.^2);
                        % Calculate the violation threshold
                        CVO = sum(max(0,Offspring.con));
                        CV  = sum(max(0,Population(List).cons),2);
                        tau = mean(CV)*sum(CV==0)/length(CV);
                        % Replace one parent with the offspring
                        replace = (d2_new<d2_old|d2_new==d2_old&d1_new<d1_old) & (CVO<tau&CV<tau|CVO==CV) | (CVO>=tau&CVO<CV);
                        Population(List(find(replace,1))) = Offspring;
                        % Update the intercepts
                        a = Intercepts(Population.objs);
                        % Update the ideal point
                        z = min(z,Offspring.obj);
                    end
                end
            end
        end
    end
end