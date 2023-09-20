classdef module_cre < handle
%module_cre - Creation module.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        GUI;                % The GUI object
        app = struct();     % All the components
        Blocks;             % All blocks
        Graph;              % Adjacent matrix
        Boxes;              % Boxes for all blocks
        Lines;              % Lines for all edges
    end
    methods(Access = ?GUI)
        %% Constructor
        function obj = module_cre(GUI)
            % The main grid
            obj.GUI = GUI;
            obj.app.maingrid = GUI.APP(3,1,uilabel(obj.GUI.app.maingrid,'Text','Under development','HorizontalAlignment','center'));
        end
    end
end