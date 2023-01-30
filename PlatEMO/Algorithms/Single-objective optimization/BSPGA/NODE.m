classdef NODE < handle
% A node in the binary space partition tree

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        node;   % Decision vector
        left;   % Left child
        right;  % Right child
        level;  % Level
    end
    methods
        function obj = NODE(varargin)
            [obj.node,obj.level] = deal(varargin{:});
        end
    end
end