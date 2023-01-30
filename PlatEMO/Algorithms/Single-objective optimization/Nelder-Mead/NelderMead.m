classdef NelderMead < ALGORITHM
% <single> <real>
% The Nelder-Mead algorithm

%------------------------------- Reference --------------------------------
% J. A. Nelder and R. Mead, A simplex method for function minimization,
% Computer Journal, 1965, 7(4): 308-313.
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
            %% Generate random solutions
            P = Problem.Initialization(1);
            P = [P,Problem.Evaluation(repmat(P.dec,Problem.D,1).*(1+eye(Problem.D)*0.05))];

            %% Optimization
            while Algorithm.NotTerminated(P)
                [~,rank] = sort(P.objs);
                P        = P(rank);
                Mean     = mean(P(1:end-1).decs,1);
                R        = Problem.Evaluation(2*Mean-P(end).dec);
                if P(1).obj <= R.obj && R.obj < P(end-1).obj
                    P(end) = R;
                    continue;
                elseif R.obj < P(1).obj
                    S = Problem.Evaluation(Mean+2*(Mean-P(end).dec));
                    if S.obj < R.obj
                        P(end) = S;
                    else
                        P(end) = R;
                    end
                    continue;
                elseif P(end-1).obj <= R.obj && R.obj < P(end).obj
                    C = Problem.Evaluation(Mean+(R.dec-Mean)/2);
                    if C.obj < R.obj
                        P(end) = C;
                        continue;
                    end
                else
                    C = Problem.Evaluation(Mean+(P(end).dec-Mean)/2);
                    if C.obj < P(end).obj
                        P(end) = C;
                        continue;
                    end
                end
                P(2:end) = Problem.Evaluation(repmat(P(1).dec,Problem.D,1)+(P(2:end).decs-repmat(P(1).dec,Problem.D,1))/2);
            end
        end
    end
end