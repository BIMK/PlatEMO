classdef SHADE < ALGORITHM
% <single> <real/integer> <large/none> <constrained/none>
% Success-history based adaptive differential evolution

%------------------------------- Reference --------------------------------
% R. Tanabe and A. Fukunaga, Success-history based parameter adaptation for
% differential evolution, Proceedings of the IEEE Congress on Evolutionary
% Computation, 2013.
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
            %% Generate random population
            Population = Problem.Initialization();
            Archive    = [];
            MCR = zeros(Problem.N,1) + 0.5;
            MF  = zeros(Problem.N,1) + 0.5;
            k   = 1;
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                % Generate offspring
                [~,rank] = sort(FitnessSingle(Population));
                Xpb = Population(rank(ceil(rand(1,Problem.N).*max(2,rand(1,Problem.N)*0.2*Problem.N)))).decs;
                Xr1 = Population(randi(end,1,Problem.N)).decs;
                P   = [Population,Archive];
                Xr2 = P(randi(end,1,Problem.N)).decs;
                CR  = randn(Problem.N,1).*sqrt(0.1) + MCR(randi(end,Problem.N,1));
                CR  = repmat(max(0,min(1,CR)),1,Problem.D);
                F   = min(1,trnd(1,Problem.N,1).*sqrt(0.1) + MF(randi(end,Problem.N,1)));
                while any(F<=0)
                    F(F<=0) = min(1,trnd(1,sum(F<=0),1).*sqrt(0.1) + MF(randi(end,sum(F<=0),1)));
                end
                F = repmat(F,1,Problem.D);
                Site         = rand(size(CR)) < CR;
                OffDec       = Population.decs;
                OffDec(Site) = OffDec(Site) + F(Site).*(Xpb(Site)-OffDec(Site)+Xr1(Site)-Xr2(Site));
                Offspring    = Problem.Evaluation(OffDec);
                % Update the population and archive
                delta   = FitnessSingle(Population) - FitnessSingle(Offspring);
                replace = delta > 0;
                Archive = [Archive,Population(replace)];
                Archive = Archive(randperm(end,min(end,Problem.N)));
                Population(replace) = Offspring(replace);
                % Update the parameters
                if any(replace)
                    w      = delta(replace)./sum(delta(replace));
                    MCR(k) = w'*CR(replace,1);
                    MF(k)  = (w'*F(replace,1).^2)./(w'*F(replace,1));
                    k      = mod(k,length(MCR)) + 1;
                end
            end
        end
    end
end