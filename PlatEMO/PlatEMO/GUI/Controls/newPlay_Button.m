classdef newPlay_Button < newButton
%newPlay_Button - The object used in newPlay.
%
%   This is the class of slipper used in newPlay, which cannot be
%   instantiated independently.
%
%   See also newPlay

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = ?newPlay)
        newvalue  = 0;  % The new value for newPlay (delayed assignment)
        presstime = 0;  % The time when this button is pressed
    end
    methods(Access = ?newPlay)
        %% Constructor
        function obj = newPlay_Button(parent,pos)
            obj@newButton(parent,pos,'| | |','FontSize',6);
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
            obj.position(1) = min(max(obj.position(1)+offset(1),1),obj.parent.labelLine.position(3)-obj.position(3)+1);
            obj.newvalue    = (obj.position(1)-1)./(obj.parent.position(3)-obj.position(3));
            if cputime-obj.presstime>=0.1 || ((obj.newvalue<=0||obj.newvalue>=obj.parent.maxvalue)&&obj.parent.value~=obj.newvalue)
                obj.parent.value = obj.newvalue;
            end
        end
    end
end