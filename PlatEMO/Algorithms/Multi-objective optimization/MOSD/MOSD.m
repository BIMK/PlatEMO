classdef MOSD < ALGORITHM
% <multi> <real> <large/none> <constrained/none>
% Multiobjective steepest descent
% step --- 0.1 --- Step size

%------------------------------- Reference --------------------------------
% X. Liu and A. C. Reynolds, A multiobjective steepest descent method with
% applications to optimal well control, Computational Geosciences, 2016,
% 20: 355-374.
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
            step = Algorithm.ParameterSet(0.1);
            
            %% Generate random population
            Population = Problem.Initialization(Problem.N);
            Archive    = UpdateArchive(Population,Problem.N);
            step0      = step;
            
            %% Optimization
            while Algorithm.NotTerminated(Archive)
                for i = 1 : Problem.N
                    gk  = FiniteDifference(Problem,Population(i)); 
                    gk1 = norm(gk(:,1));
                    gk2 = norm(gk(:,2));
                    if gk1 <= gk2
                        gs = gk(:,1);
                        gl = gk(:,2);
                    else
                        gs = gk(:,2);
                        gl = gk(:,1);
                    end
                    if gs'*(gs-gl) <= 0
                        d = -gs/norm(gs);
                    else
                        D = (((norm(gl).^2-gs'*gl)/(norm(gs).^2-gs'*gl))*(-gs))-gl;
                        d = D/norm(D);
                    end
                    Offspringdec = Population(i).dec+step*d';
                    Offspring    = Problem.Evaluation(Offspringdec);
                    Archive      = UpdateArchive([Archive,Offspring],Problem.N);
                    if ~any(Offspring.obj<Population(i).obj)
                        Population(i) = Archive(randi(end));
                        step = step/2;
                    else
                        Population(i) = Offspring;
                        step = min([2*step,step0]);
                    end    
                end              
            end
        end
    end
end

function df = FiniteDifference(Problem,X)
    if any(X.con>0)
        df = Problem.CalConGrad(X.dec)';
    else
        df = Problem.CalObjGrad(X.dec)';
    end
end

function P = UpdateArchive(P,N)
    P = P(NDSort(P.objs,P.cons,1)==1);
    if length(P) > N
        Choose = true(1,length(P));
        Dis    = pdist2(P.objs,P.objs);
        Dis(logical(eye(length(Dis)))) = inf;
        while sum(Choose) > N
            Remain   = find(Choose);
            Temp     = sort(Dis(Remain,Remain),2);
            [~,Rank] = sortrows(Temp);
            Choose(Remain(Rank(1))) = false;
        end
        P = P(Choose);
    end
end