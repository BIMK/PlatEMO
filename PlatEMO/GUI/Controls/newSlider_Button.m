classdef newSlider_Button < newPanel
%newSlider_Button - The object used in newSlider.
%
%   This is the class of slipper used in newSlider, which cannot be
%   instantiated independently.
%
%   See also newSlider

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = ?newSlider)
        newvalue  = 0;  % The new value for newSlider (delayed assignment)
        presstime = 0;  % The time when this button is pressed
    end
    methods(Access = ?newSlider)
        %% Constructor
        function obj = newSlider_Button(parent,pos)
            obj = obj@newPanel(parent,pos,[]);
            % Update the button style by obj.moved, obj.pressed, obj.value
            % or obj.state
            obj.newListener(obj,{'moved','pressed','state'},'PostSet',@obj.updateStyle);
        end
    end
    methods(Access = ?newGUI)
        function mouseDown(obj)
            obj.presstime = cputime;
        end
        function mouseUp(obj)
            if obj.parent.value ~= obj.newvalue
                obj.parent.value = obj.newvalue;
            end
        end
        function moveOut(obj)
            if obj.parent.value ~= obj.newvalue
                obj.parent.value = obj.newvalue;
            end
        end
        function mouseDrag(obj,offset)
            obj.position(2) = min(max(obj.position(2)+offset(2),16),obj.parent.position(4)-54);
            obj.newvalue    = 1 - (obj.position(2)-16)./(obj.parent.position(4)-70);
            if cputime-obj.presstime>=0.1 || ((obj.newvalue<=0||obj.newvalue>=1)&&obj.parent.value~=obj.newvalue)
                obj.parent.value = obj.newvalue;
            end
        end
    end
    methods(Access = protected)
        %% Update the button style
        function updateStyle(obj,hObject,eventdata)
            if obj.state
                if obj.pressed
                    obj.handle.BackgroundColor = (obj.parent.parent.color+[.6 .6 1])/2;
                elseif obj.moved
                    obj.handle.BackgroundColor = (obj.parent.parent.color+[.8 .8 1])/2;
                else
                    obj.handle.BackgroundColor = obj.parent.parent.color;
                end
            else
                obj.handle.BackgroundColor = obj.parent.parent.color;
            end
        end
    end
end