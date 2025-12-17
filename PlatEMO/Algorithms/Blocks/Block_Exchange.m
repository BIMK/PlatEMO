classdef Block_Exchange < BLOCK
% Exchange of parents
% nParents --- 2 --- Number of parents generating one offspring

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        nParents;   % <hyperparameter> Number of parents for generating an offspring
        Fitness;    % <parameter> Probability of selecting each parent
    end
    methods
        %% Default settings of the block
        function obj = Block_Exchange(nParents)
            obj.nParents = nParents;	% Hyperparameter
            obj.lower    = zeros(1,nParents) + 1e-20;   % Lower bounds of parameters
            obj.upper    = ones(1,nParents);            % Upper bounds of parameters
            % Randomly set the parameters
            obj.parameter = unifrnd(obj.lower,obj.upper);
            obj.ParameterAssign();
        end
        %% Assign parameters to variables
        function ParameterAssign(obj)
            obj.Fitness = cumsum(obj.parameter);
            obj.Fitness = obj.Fitness./max(obj.Fitness);
        end
        %% Main procedure of the block
        function Main(obj,Problem,Precursors,Ratio)
            ParentDec  = obj.Gather(Problem,Precursors,Ratio,2,obj.nParents);
            selected   = arrayfun(@(S)find(rand<=obj.Fitness,1),zeros(size(ParentDec,1)/obj.nParents,size(ParentDec,2)));
            selected   = selected + repmat((0:size(selected,1)-1)'*obj.nParents,1,size(ParentDec,2));
            selected   = selected + repmat((0:size(selected,2)-1)*size(ParentDec,1),size(selected,1),1);
            obj.output = ParentDec(selected);
        end
    end
end