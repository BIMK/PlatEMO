classdef TELSO < ALGORITHM
% <2024> <multi> <real/binary> <large/none> <constrained/none> <sparse>
% Two-layer encoding learning swarm optimizer

%------------------------------- Reference --------------------------------
% S. Qi, R. Wang, T. Zhang, X. Yang, R. Sun, and L. Wang. A two-layer
% encoding learning swarm optimizer based on frequent itemsets for sparse
% large-scale multi-objective optimization. IEEE/CAA Journal of Automatica
% Sinica, 2024, 11(6): 1342-1357.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Sheng Qi (email: 2745679162@qq.com)

    methods
        function main(Algorithm,Problem)
            %% Initialization reference vector
            [V,Problem.N] = UniformPoint(Problem.N,Problem.M);
            
            %% Generate random population
            % Calculate the fitness of each decision variable
            TDec    = [];
            TMask   = [];
            TempPop = [];
            DF      = zeros(1,Problem.D);	% The fitness of decision variables
            REAL    = all(Problem.encoding==1);
            for i = 1 : 1+4*REAL
                if REAL
                    Dec = unifrnd(repmat(Problem.lower,Problem.D,1),repmat(Problem.upper,Problem.D,1));
                else
                    Dec = ones(Problem.D,Problem.D);
                end
                Mask       = eye(Problem.D);
                Population = Problem.Evaluation(Dec.*Mask);
                TDec       = [TDec;Dec];
                TMask      = [TMask;Mask];
                TempPop    = [TempPop,Population];
                DF         = DF + NDSort([Population.objs,Population.cons],inf);
            end
            
            %% Generate initial population
            if REAL
                Dec = unifrnd(repmat(Problem.lower,Problem.N,1),repmat(Problem.upper,Problem.N,1));
            else
                Dec = ones(Problem.N,Problem.D);
            end
            Mask = zeros(Problem.N,Problem.D);
            for i = 1 : Problem.N
                Mask(i,TournamentSelection(2,ceil(rand*Problem.D),DF)) = 1;
            end
            Population = Problem.Evaluation(Dec.*Mask);
            [Population,Mask] = Initial_EnvironmentalSelection([Population,TempPop],[Dec;TDec],[Mask;TMask],Problem.N);
            [Population,Mask] = EnvironmentalSelection(Population,V,(Problem.FE/Problem.maxFE)^2,Mask);
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                Fitness    = CalFitness(Population.objs);
                [~,index]  = sort(Fitness);
                Population = Population(index);
                Mask       = Mask(index,:);
                [Offspring,Mask]  = Operator(Population,Mask,Problem);
                [Population,Mask] = EnvironmentalSelection([Population,Offspring],V,(Problem.FE/Problem.maxFE)^2,Mask);
            end
        end
    end
end

function Fitness = CalFitness(PopObj)
% Calculate the fitness by shift-based density
    N      = size(PopObj,1);
    fmax   = max(PopObj,[],1);
    fmin   = min(PopObj,[],1);
    PopObj = (PopObj-repmat(fmin,N,1))./repmat(fmax-fmin,N,1);
    Dis    = inf(N);
    for i = 1 : N
        SPopObj = max(PopObj,repmat(PopObj(i,:),N,1));
        for j = [1:i-1,i+1:N]
            Dis(i,j) = norm(PopObj(i,:)-SPopObj(j,:));
        end
    end
    Fitness = min(Dis,[],2);
end