classdef newLabelButton < newGUI
%newLabelButton - Label button.
%
%   h = newLabelButton(Parent,Pos,Str,...) creates a newLabel-based button
%   which has a parent of Parent, a position of Pos and a text of Str.
%
%   h.value denotes the times of clicking on the button. If the number of
%   times is even, then h.value = false, otherwise h.value = true. This
%   property is meaningful when h.choosed = true.
%
%   h.choosed = true denotes that it is a toggle button, otherwise it is a
%   general button.
%
%   The newLabelButton is similar to newButton besides a frame, and the
%   former has a simpler framework than the latter.
%
%   Example:
%       newLabelButton(h,[10 10 100 100,1 0 1 0],'test')
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

    properties(SetObservable)
        value   = false;    % Whether the button is being chosen
        choosed = false;    % Whether the button can be chosen
    end
    properties(Access = protected)
        foregroundcolor;    % Foreground color of the object
    end
    methods
        %% Constructor
        function obj = newLabelButton(parent,pos,str,varargin)
            handle = uicontrol('FontName','Microsoft YaHei','FontSize',11,'BackgroundColor',parent.color,...
                               'Parent',parent.handle,'Style','text','String',str,'Enable','inactive','Units','pixels','Position',pos(1:4));
            obj@newGUI(3,parent,handle,pos(5:8),varargin{:});
            obj.foregroundcolor = obj.handle.ForegroundColor;
            % Update the button style by obj.moved, obj.pressed, obj.value
            % or obj.state
            obj.newListener(obj,{'moved','pressed','value','state'},'PostSet',@obj.updateStyle);
        end
    end
    methods(Access = ?newGUI)
        function mouseUp(obj)
            obj.value = ~obj.value;
            obj.callback(obj,[]);
        end
    end
    methods(Access = protected)
        %% Update the button style
        function updateStyle(obj,hObject,eventdata)
            if obj.state
            	obj.handle.ForegroundColor = obj.foregroundcolor;
                if obj.value && obj.choosed || obj.pressed
                    obj.handle.BackgroundColor = (obj.parent.color+[.6 .6 1])/2;
                elseif obj.moved
                    obj.handle.BackgroundColor = (obj.parent.color+[.8 .8 1])/2;
                else
                    obj.handle.BackgroundColor = obj.parent.color;
                end
            else
                obj.handle.ForegroundColor = (obj.foregroundcolor+obj.parent.color)/2;
                if obj.value && obj.choosed
                    obj.handle.BackgroundColor = [.9 .9 .9];
                else
                    obj.handle.BackgroundColor = obj.parent.color;
                end
            end
        end
    end
end