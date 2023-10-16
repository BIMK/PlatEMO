classdef MOCGDE < ALGORITHM
% <multi/many> <real> <large/none> <constrained/none>
% Multi-objective conjugate gradient and differential evolution algorithm
% NP --- 10 --- Small population size

%------------------------------- Reference --------------------------------
% Y. Tian, H. Chen, H. Ma, X. Zhang, K. C. Tan, and Y. Jin, Integrating
% conjugate gradients into evolutionary algorithms for large-scale
% continuous multi-objective optimization, IEEE/CAA Journal of Automatica
% Sinica, 2022, 9(10): 1801-1817.
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
            %% Generate the weight vectors
            NP       = Algorithm.ParameterSet(10);
            [W,subN] = UniformPoint(NP,Problem.M);
            
            %% Generate random population
            Population = Problem.Initialization(subN);
            Archive    = Population;
            
            %% Optimization
            K  = zeros(1,subN);
            g0 = cell(1,subN);
            d0 = cell(1,subN);
            while Algorithm.NotTerminated(Archive)
                K = mod(K,Problem.D) + 1;
                OffPop = [];
                for i = 1 : subN
                    [gk,site] = FiniteDifference(Problem,Population(i),W(i,:));
                    if K(i) == 1
                        dk = -gk;
                    else
                        beta = (gk'*gk)/(g0{i}'*g0{i});
                        dk   = -gk + beta*d0{i};
                        if gk'*dk >= 0
                            dk = -gk;
                        end
                    end
                    success = false;
                    for step = 0 : 9
                        mu        = rand(1,Problem.D) < 1/sum(site);
                        OffDec    = Population(i).dec + ~site.*0.5^step.*dk' + mu.*site*0.5^step.*(Archive(randi(end)).dec-Archive(randi(end)).dec);
                        Offspring = Problem.Evaluation(OffDec);
                        OffPop    = [OffPop,Offspring];
                        if sum(max(Offspring.con,0))<sum(max(Population(i).con,0)) || sum(max(Offspring.con,0))==sum(max(Population(i).con,0)) && all(Offspring.obj<Population(i).obj)
                            success = true;
                            break;
                        end
                    end
                    if success
                        Population(i) = Offspring;
                        g0{i} = gk;
                        d0{i} = dk;
                    else
                        Population(i) = Archive(randi(end));
                        K(i) = 0;
                    end
                end
                Archive = UpdateArchive([Archive,OffPop],Problem.N);
            end
        end
    end
end

function [df,site] = FiniteDifference(Problem,X,W)
    if any(X.con>0)
        df   = Problem.CalConGrad(X.dec)';
        site = false(1,length(X.dec));
        df   = sum(df,2);
    else
        df   = Problem.CalObjGrad(X.dec)';
        site = (any(df<0,2)&any(df>0,2))';
        df   = df*W';
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