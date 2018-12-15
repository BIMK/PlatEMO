classdef newPopmenu2_Item < newLabelButton
%newPopmenu2_Item - The object used in newPopmenu2.
%
%   This is the class of item used in newPopmenu2, which cannot be
%   instantiated independently.
%
%   See also newPopmenu2

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods(Access = ?newPopmenu2)
        %% Constructor
        function obj = newPopmenu2_Item(parent,pos,str,varargin)
            obj@newLabelButton(parent,pos,str,varargin{:});
        end
    end
    methods(Access = ?newGUI)
        function mouseUp(obj)
            obj.parent.parent.index = find([obj.parent.parent.items.handle]==obj.handle,1);
        end
    end
    methods(Access = protected)
        %% Update the button style
        function updateStyle(obj,hObject,eventdata)
            if obj.state
                if obj.pressed
                    obj.handle.BackgroundColor = (obj.parent.color+[.6 .6 1])/2;
                elseif obj.moved
                    obj.handle.BackgroundColor = (obj.parent.color+[.8 .8 1])/2;
                else
                    obj.handle.BackgroundColor = obj.parent.color;
                end
            else
                obj.handle.BackgroundColor = obj.parent.color;
            end
            if ~isempty(obj.parent.parent.delButton)
                if obj.moved
                    obj.parent.parent.delButton.visible     = true;
                    obj.parent.parent.delButton.position(2) = obj.position(2);
                    obj.parent.parent.delButton.callback    = @(~,~)obj.parent.parent.del(find([obj.parent.parent.items.handle]==obj.handle,1));
                elseif obj.parent.parent.delButton.position(2) == obj.position(2)
                    obj.parent.parent.delButton.visible = false;
                end
            end
        end
    end
end