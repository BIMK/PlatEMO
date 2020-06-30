classdef newSlider < newPanel
%newSlider - Slider.
%
%   h = newSlider(Parent,Pos,...) creates a vertical slider which has a
%   parent of Parent and a position of Pos.
%
%   h.value denotes the current value of the slider, which has a minimum
%   value of 0 and a maximum value of 1. h.callback will be invoked when
%   h.value is set.
%
%   Example:
%       newSlider(f,[10 10 100 100,1 0 1 0])
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
        value = 0;                  % The current value of the bar
    end
    properties(SetAccess = protected)
        button;                     % The three buttons on the bar
    end
    methods
        %% Constructor
        function obj = newSlider(parent,pos,varargin)
            pos(3) = 16;
            obj@newPanel(parent,pos,parent.color,varargin{:});
            % Create a background
            newLabel(obj,[1 1 pos(3:4),1 1 1 1],'','BackgroundColor',([.9 .9 .9]+parent.color)/2,'callback',@obj.cb_changeValue);
            % Create the three buttons on the slider
            obj.button    = newButton(obj,[1 pos(4)-14 pos(3) 15,1 1 0 1],'¡ø','ForegroundColor',[0.5 0.5 0.5],'FontSize',8,'callback',@(~,~)obj.set('value',obj.value-0.1));
            obj.button(2) = newButton(obj,[1 1 pos(3) 15,1 1 1 0],'¨‹','ForegroundColor',[0.5 0.5 0.5],'FontSize',8,'callback',@(~,~)obj.set('value',obj.value+0.1));
            obj.button(3) = newSlider_Button(obj,[1 pos(4)-54 pos(3) 40,1 1 1 0]);
            % Update the position of obj.button(3) by obj.handle.position
            obj.newListener(obj.handle,'SizeChanged',@obj.cb_updateButton);
        end
        %% Update the position of obj.button(3) by obj.value
        function set.value(obj,value)
            obj.value = min(max(value,0),1);
            obj.button(3).position(2) = round(16+(obj.position(4)-70)*(1-obj.value));
            obj.button(3).newvalue    = obj.value;
            obj.button(3).presstime   = cputime;
            obj.callback(obj,obj.value);
        end
    end
    methods(Access = ?newGUI)
        function mouseScroll(obj,moveDown)
            if moveDown
                obj.value = obj.value + 0.1;
            else
                obj.value = obj.value - 0.1;
            end
        end
    end
    methods(Access = protected)
        %% Change obj.value when click on the bar
        function cb_changeValue(obj,hObject,eventdata)
            obj.value = 1 - (obj.figure.handle.CurrentPoint(2)-obj.absPos(2)-16)./(obj.position(4)-30);
        end
        %% Change obj.value when its size is changed
        function cb_updateButton(obj,hObject,eventdata)
            obj.button(3).position(2) = 16 + (obj.position(4)-70)*(1-obj.value);
        end
    end
end