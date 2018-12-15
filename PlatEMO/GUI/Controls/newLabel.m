classdef newLabel < newGUI
%newLabel - Label.
%
%   h = newLabel(Parent,Pos,Str,...) creates a label which has a parent of
%   Parent, a position of Pos and a text of Str.
%
%   Example:
%       newLabel(f,[10 10 100 100,1 0 1 0],'test')
%
%   See also newGUI

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

	properties(SetAccess = protected)
        foregroundcolor;    % Foreground color of the object
        backgroundcolor;    % Background color of the object
    end
    methods
        %% Constructor
        function obj = newLabel(parent,pos,str,varargin)
            handle = uicontrol('FontName','Microsoft YaHei','FontSize',11,'BackgroundColor',parent.color,...
                               'Parent',parent.handle,'Style','text','String',str,'Enable','inactive','Units','pixels','Position',pos(1:4));
            obj@newGUI(3,parent,handle,pos(5:8),varargin{:});
            obj.foregroundcolor = obj.handle.ForegroundColor;
            obj.backgroundcolor = obj.handle.BackgroundColor;
            % Update the label style by obj.state
            obj.newListener(obj,'state','PostSet',@obj.updateStyle);
        end
    end
    methods(Access = protected)
        %% Update the label style
        function updateStyle(obj,hObject,eventdata)
            if obj.state
                obj.handle.ForegroundColor = obj.foregroundcolor;
                obj.handle.BackgroundColor = obj.backgroundcolor;
            else
                obj.handle.ForegroundColor = (obj.foregroundcolor+obj.parent.color)/2;
                obj.handle.BackgroundColor = (obj.backgroundcolor+obj.parent.color)/2;
            end
        end
    end
end