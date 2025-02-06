classdef GEA < ALGORITHM
% Graph evolutionary algorithm

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [Blocks,Graph] = Algorithm.ParameterSet({},{});
            
            %% Generate random population
            isPop = arrayfun(@(s)isa(s,'Block_Population'),Blocks(:)');
            Blocks(isPop).Initialization(Problem.Initialization());

            %% Optimization
            while Algorithm.NotTerminated(Blocks(1).output)
                activated = false(1,length(Blocks));
                while ~all(activated(isPop))
                    for i = find(~activated)
                        if all(activated(logical(Graph(:,i)))|isPop(logical(Graph(:,i))))
                        	Blocks(i).Main(Problem,Blocks(logical(Graph(:,i))),Graph(:,i));
                            activated(i) = true;
                        end
                    end
                end
            end
        end
    end
end