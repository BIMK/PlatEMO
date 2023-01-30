classdef MOEGS < ALGORITHM
% <multi> <real> <large/none>
% Multi-objective evolutionary gradient search
% r --- 500 --- Number of test candidates

%------------------------------- Reference --------------------------------
% C. K. Goh, Y. S. Ong, K. C. Tan, and E. J. Teoh, An investigation on
% evolutionary gradient search for multi-objective optimization,
% Proceedings of the IEEE Congress on Evolutionary Computation, 2008.
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
            r = Algorithm.ParameterSet(500);
            
            %% Generate random population
            Population = Problem.Initialization(Problem.N);
            Archive    = UpdateArchive(Population,Problem.N);
            thao       = 0.1;
            
            %% Optimization
            while Algorithm.NotTerminated(Archive)
                for i =1 : Problem.N
                    parent    = Population(i);
                    z         = normrnd(0,thao,[Problem.D,r])';
                    t         = repmat(parent.dec,r,1)+z;
                    Offspring = Problem.Evaluation(t);              
                    Archive   = UpdateArchive([Archive,Offspring],Problem.N);
                    if any(parent.con>0)   
                         gk = sum((Offspring.cons-parent.con),2);
                    elseif size(Offspring.objs-parent.obj,2) == 2 
                        w = rand(r,1);
                        w = [w,1-w];
                        gk = sum(w.*(Offspring.objs-parent.obj),2);  
                    elseif size(Offspring.objs-parent.obj,2) == 3
                        w  = rand(r,1);
                        w  = [w/3,(1-w/3)/2,1-w/3-(1-w/3)/2];
                        gk = sum(w.*(Offspring.objs-parent.obj),2);                    
                    end
                    c = parent.dec - thao*gk'*z/norm(gk'*z);
                    Offspring = Problem.Evaluation(c);
                    if sum(max(Offspring.con,0))<sum(max(parent.con,0)) || all(Offspring.obj<parent.obj)
                        Population(i) = Offspring;
                        thao = thao*1.8;
                    elseif sum(max(parent.con,0))<sum(max(Offspring.con,0)) || all(parent.obj<Offspring.obj)
                        Population(i) = Archive(randi(end));
                        thao = thao/1.8;
                    else
                        if rand()>0.5
                            Population(i) = Offspring;
                        else
                            Population(i) = Archive(randi(end));
                        end
                    end
                    Archive = UpdateArchive([Archive,Offspring],Problem.N);
                    if thao > 1.2
                        thao = thao/1.2;
                    elseif thao < 1e-4
                        thao = thao*1.2;
                    end
                end
            end
        end
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