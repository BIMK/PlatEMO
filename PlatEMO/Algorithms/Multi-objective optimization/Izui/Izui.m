classdef Izui < ALGORITHM
% <multi/many> <real> <large/none> <constrained/none>
% An aggregative gradient based multi-objective optimizer proposed by Izui et al.

%------------------------------- Reference --------------------------------
% K. Izui, T. Yamada, S. Nishiwaki, and K. Tanaka, Multiobjective
% optimization using an aggregative gradient-based method, Structural and
% Multidisciplinary Optimization, 2015, 51: 173-182.
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
            %% Generate random population
            Population = Problem.Initialization(Problem.N);
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                for i = 1 : Problem.N
                    w             = Calw(Population,Population(i));
                    PopObj_g      = Calg(Problem,Population(i));
                    OffspringDec  = Calf(w,PopObj_g,Problem.lower,Problem.upper);
                    Offspring     = Problem.Evaluation(OffspringDec');
                    Population(i) = Offspring;
                end
            end
        end
    end
end

function df = Calg(Problem,X)
    if any(X.con>0)
        df = Problem.CalConGrad(X.dec)';
    else
        df = Problem.CalObjGrad(X.dec)';
    end
end

function w = Calw(Population,X)
    f   = X.objs;
    A   = Population.objs;
    b   = ones(1,length(Population));
    Aeq = [];
    beq = [];
    lb  = zeros(1,length(f));
    op  = optimoptions('linprog','Display','none');
    w   = linprog(f,-A,-b,Aeq,beq,lb,[],op);
end

function OffspringDec = Calf(w,PopObj_g,lower,upper)
    f   = (PopObj_g*w)';
    A   = [];
    b   = [];
    Aeq = [];
    beq = [];
    lb  = 0.05*lower;
    ub  = 1.05*upper;
    op  = optimoptions('linprog','Display','none');
    OffspringDec = linprog(f,A,b,Aeq,beq,lb,ub,op);
end