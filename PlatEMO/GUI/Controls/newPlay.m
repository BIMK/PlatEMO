classdef newPlay < newPanel
%newPlay - Play bar.
%
%   h = newPlay(Parent,Pos,...) creates a horizontal play bar which has a
%   parent of Parent and a position of Pos.
%
%   h.value denotes the current value of the bar and h.maxvalue denotes the
%   maximum value of h.value can be set, both of which have a minimum value
%   of 0 and a maximum value of 1. h.callback will be invoked when h.value
%   is set.
%
%   Example:
%       newPlay(f,[10 10 100 100,1 0 1 0])
%
%   See also newGUI, newPanel

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetObservable)
        value    = 0;               % The current value of the bar
        maxvalue = 0;               % The maximal value of the bar
    end
    properties(SetAccess = protected)
        labelLine;                  % The process bar of the object
        button;                     % The slipper button on the object
    end
    methods
        %% Constructor
        function obj = newPlay(parent,pos,varargin)
            obj@newPanel(parent,pos,parent.color,varargin{:});
            % Create a background
            x = newLabel(obj,[1 pos(4)*0.36 pos(3) pos(4)*0.33,1 1 0 0],'','BackgroundColor',[.95 .95 .95],'callback',@obj.mouseUp);
            % Create the process bar
            obj.labelLine = newLabel(obj,[1 pos(4)*0.36 30 pos(4)*0.33,1 0 0 0],'','BackgroundColor',[0.078 0.169 0.549],'callback',@obj.mouseUp);
            % Create the slipper
            obj.button = newPlay_Button(obj,[1 pos(4)*0.15 30 pos(4)*0.7,1 0 1 1]);
            % Update the position of obj.button by obj.handle.position
            obj.newListener(obj.handle,'SizeChanged',@obj.cb_updateButton);
        end
        %% The set method of obj.value
        function set.value(obj,value)
            obj.value = min(max(value,0),obj.maxvalue);
            obj.button.position(1) = round(1+(obj.position(3)-obj.button.position(3))*obj.value);
            obj.button.newvalue    = obj.value;
            obj.button.presstime   = cputime;
            obj.callback(obj,obj.value);
        end
        %% The set method of obj.maxvalue
        function set.maxvalue(obj,value)
            obj.maxvalue = min(max(value,0),1);
            obj.labelLine.position(3) = obj.button.position(3) + (obj.position(3)-obj.button.position(3))*obj.maxvalue;
            if obj.value > obj.maxvalue
                obj.value = obj.maxvalue;
            end
        end
    end
    methods(Access = ?newGUI)
        function mouseScroll(obj,moveDown)
            if moveDown
                obj.value = obj.value - 0.05;
            else
                obj.value = obj.value + 0.05;
            end
        end
        function mouseUp(obj,hObject,eventdata)
            obj.value = (obj.figure.handle.CurrentPoint(1)-obj.absPos(1))./(obj.position(3)-obj.button.position(3));
        end
    end
    methods(Access = protected)
        %% Update the position of obj.button by obj.handle.position
        function cb_updateButton(obj,hObject,eventdata)
            obj.button.position(1)    = 1 + (obj.position(3)-obj.button.position(3))*obj.value;
            obj.labelLine.position(3) = obj.button.position(3) + (obj.position(3)-obj.button.position(3))*obj.maxvalue;
        end
    end
end