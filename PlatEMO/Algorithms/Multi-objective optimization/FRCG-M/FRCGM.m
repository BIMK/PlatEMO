classdef FRCGM < ALGORITHM
% <multi/many> <real> <large/none> <constrained/none>
% Fletcher-Reeves conjugate gradient (for multi-objective optimization)
% beta  --- 0.6 --- A parameter within [0,1] for line search
% sigma --- 0.4 --- A parameter within [0 0.5] for line search

%------------------------------- Reference --------------------------------
% R. Fletcher and C. M. Reeves, Function minimization by conjugate
% gradients, The Computer Journal, 1964, 7(2): 149-154.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [beta,sigma] = Algorithm.ParameterSet(0.6,0.4);
            
            %% Generate random population
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            Population    = Problem.Initialization();

            %% Optimization
            k  = 0;
            g0 = cell(1,Problem.N);
            d0 = cell(1,Problem.N);
            while Algorithm.NotTerminated(Population)
                for i = 1 : Problem.N                   
                    gk    = FiniteDifference(Problem,Population(i),W(i,:));
                    itern = k - (Problem.D+1)*floor(k/(Problem.D+1)) + 1;
                    if itern <= 1
                        dk = -gk;
                    else
                        betak = (gk'*gk)/(g0{i}'*g0{i});
                        dk    = -gk + betak*d0{i};
                        gd    = gk'*dk;
                        if gd >= 0
                            dk = -gk;
                        end
                    end
                    for m = 0 : 20
                        X = Problem.Evaluation(Population(i).dec+beta^m*dk');
                        if X.obj*W(i,:)' <= Population(i).obj*W(i,:)' + sigma*beta^m*gk'*dk
                            break;
                        end
                    end
                    Population(i) = X;
                    g0{i} = gk;
                    d0{i} = dk;
                end
                k  = k + 1;
            end
        end
    end
end

function df = FiniteDifference(Problem,X,W)
    if any(X.con>0)
        df = Problem.CalConGrad(X.dec)';
        df = sum(df,2);
    else
        df = Problem.CalObjGrad(X.dec)';
        df = df*W';
    end
end