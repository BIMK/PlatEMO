classdef newTip < newButtonSpecial
%newTip - Tip button.
%
%   h = newTip(Parent,Pos,O,...) creates a tip button which has a parent of
%   Parent and a position of Pos, and it is associated with control O. That
%   is, h always has the same state to O, and locates to the right of O.
%
%   Example:
%       newTip(f,[10 10 100 100,1 0 1 0],Obj)
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
        targetObj;      % Target object to be associated with
    end
    methods
        %% Constructor
        function obj = newTip(parent,pos,targetObj,varargin)
            % Create the button
            obj@newButtonSpecial(parent,[1 1 pos(3) targetObj.handle.Position(4),1 0 floor(targetObj.sizeLock(2)/10) mod(targetObj.sizeLock(2),10)],2,[],varargin{:});
            % Update the button style according to targetObj
            obj.targetObj = targetObj;
            obj.newListener(obj.targetObj,{'moved','pressed','value','choosed'},'PostSet',@obj.updateStyle);
            % Update the position according to targetObj
            obj.newListener(obj.targetObj,'position','PostSet',@obj.updatePosition);
            obj.updateCData();
            obj.updatePosition();
        end
    end
    methods(Access = protected)
        %% Update the button style according to obj.targetObj
        function updateStyle(obj,hObject,eventdata)
            value   = obj.targetObj.value;
            choosed = obj.targetObj.choosed;
            pressed = obj.pressed | obj.targetObj.pressed;
            moved   = obj.moved | obj.targetObj.moved;
            if obj.state
                if value && choosed || pressed
                    obj.handle.CData = obj.CDataPress;
                elseif moved
                    obj.handle.CData = obj.CDataMove;
                else
                	obj.handle.CData = obj.CDataNormal;
                end
            else
                if value && choosed
                    obj.handle.CData = obj.CDataUnableP;
                else
                    obj.handle.CData = obj.CDataUnable;
                end
            end
        end
        %% Update the position according to obj.targetObj
        function updatePosition(obj,hObject,eventdata)
            obj.position(1:2) = [obj.targetObj.position(1)+obj.targetObj.position(3),obj.targetObj.position(2)];
        end
    end
end