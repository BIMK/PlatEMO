classdef POCEA < ALGORITHM
% <multi> <real/integer> <large/none> <constrained>
% Paired offspring generation based constrained evolutionary algorithm

%------------------------------- Reference --------------------------------
% C. He, R. Cheng, Y. Tian, X. Zhang, K. C. Tan, and Y. Jin, Paired
% offspring generation for constrained large-scale multiobjective
% optimization, IEEE Transactions on Evolutionary Computation, 2021, 25(3):
% 448-462.
%--------------------------------------------------------------------------

% Copyright (c) 2020-2021 Cheng He

    methods
        function main(Algorithm,Problem)
            %% Parameter settings
            k = Algorithm.ParameterSet(5);
            Population     = Problem.Initialization();
            [V0,Problem.N] = UniformPoint(Problem.N,Problem.M);
            [Vs0,L]        = UniformPoint(floor(Problem.N/k),Problem.M);
            [V,Vs]         = deal(V0,Vs0);
            %% Optimization
            while Algorithm.NotTerminated(Population)
                [INDEX,THETA,DIS] = Association(Population,Vs,k);
                CV = sum(max(Population.cons,0),2);
                rf = sum(CV<1e-6)/length(Population);
                Offspring = [];

                %% Paired offspring generation
                for i = 1 : L
                    % Subpopulation construction
                    SubPop = Population(INDEX(:,i));
                    theta  = THETA(INDEX(:,i));
                    if mean(theta) >= pi/L/2
                        [~,index] = sort(DIS);
                        selected  = index(1:min(k,length(index)));
                        SubPop    = [SubPop,Population(selected)];
                        epsilon   = max(CV([INDEX(:,i);selected]));
                    else
                        epsilon   = min(CV(INDEX(:,i)))*(1-rf)+mean(CV(INDEX(:,i)))*rf;
                    end
                    % Offspring generation
                    if length(SubPop) < 2
                        rank = [1, 1];
                    else
                        [~,rank]= sort(rand(k,length(SubPop)),2);
                    end
                    [winner,loser] = CHP(SubPop(rank(1)),SubPop(rank(2)),epsilon);
                    Offspring      = [Offspring,Operator(Problem,loser,winner)];
                end 
                Population = RVEASelection([Population,Offspring],V,Problem.N,(Problem.FE/Problem.maxFE)^2);
                if ~mod(ceil(Problem.FE/Problem.N),ceil(0.1*Problem.maxFE/Problem.N))
                    [V,Vs] = ReferenceVectorAdaptation(Population.objs,V0,Vs0);
                end
            end
        end
    end
end

function [INDEX,THETA,DIS] = Association(Population,V,k)
% Associate k candidate solutions to each reference vector

    PopObj    = Population.objs - repmat(min(Population.objs,[],1),length(Population),1);
	DIS	      = sum(PopObj.^2,2);
    THETA     = acos(1 - pdist2(PopObj,V,'cosine'));
    [~,index] = sort(THETA,1);
    INDEX     = index(1:min(k,end),:);
end