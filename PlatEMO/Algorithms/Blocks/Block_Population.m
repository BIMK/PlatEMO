classdef Block_Population < BLOCK
% A population

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        %% Main procedure of the block
        function Main(obj,Problem,Precursors,Ratio)
            obj.output = obj.Gather(Problem,Precursors,Ratio,1,1);
        end
        %% Initialize the solutions
        function Initialization(obj,Population)
            [obj.output] = deal(Population);
        end
    end
end