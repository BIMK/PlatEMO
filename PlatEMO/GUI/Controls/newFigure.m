classdef newFigure < newGUI
%newFigure - Figure window.
%
%   h = newFigure(Pos,Str) creates a figure which has a position of Pos and
%   a title of Str.
%
%   h = newFigure(Pos,Str,F) creates a modal subfigure which has a position
%   specified by the location of the mouse in figure F.
%
%	To improve the efficiency, there is no delayed strategy in adding or
%	deleting offsprings, so DO NOT create or delete object in the
%	movecallback of any controls, since more than one movecallback may be
%	executed in one WindowButtonMotionFcn callback of the figure.
%
%   Example:
%       newFigure([1 1 100 100],'new figure')
%       newFigure([1 1 100 100],'subfigure',f)
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

	properties
        busy = false;    	% Whether the figure is busy
    end
    properties(SetAccess = ?newGUI)
        offsprings;         % All the objects in the figure
        offState;           % State of each offspring
        offPos;             % Position of each offspring
        offMoved;           % The value of moved of each offspring
        offPressed;         % The value of pressed of each offspring
        mouseLoc;       	% The last location of the mouse
    end
    methods
        %% Constructor
        function obj = newFigure(pos,title,mainFigure)
            if nargin > 2
                % The position of subfigure
                pos(1)   = mainFigure.handle.CurrentPoint(1) + 10;
                pos(2)   = mainFigure.handle.CurrentPoint(2) - (pos(4)+10);
                pos(1:2) = pos(1:2) + mainFigure.position(1:2);
            end
            handle = figure('MenuBar','none','ToolBar','none','NumberTitle','off','Color','white','Name',title,'Position',pos(1:4));
            if nargin > 2
                % Subfigure
                handle.WindowStyle = 'modal';
                handle.Resize = 'off';
            else
                % Main figure
                movegui(handle,'center');
            end
            obj@newGUI(1,struct('state',true),handle,[1 0 1 0]);
            obj.busy     = true;
            obj.mouseLoc = obj.handle.CurrentPoint;
            uicontrol('Position',[-1 -1 0 0]);
        end
        %% Get the background color of the figure
        function value = color(obj)
            value = obj.handle.Color;
        end
        %% The set method of obj.busy
        function set.busy(obj,value)
            obj.busy = value;
            drawnow();
            if obj.busy
                set(obj.handle,'WindowButtonMotionFcn',[],'WindowScrollWheelFcn',[],'WindowButtonDownFcn',[],...
                    'WindowButtonUpFcn',[],'SizeChangedFcn',@obj.WindowResize,'CloseRequestFcn',@obj.WindowClose);
                obj.handle.Pointer = 'watch';
            else
                set(obj.handle,'WindowButtonMotionFcn',@obj.WindowButtonMotion,'WindowScrollWheelFcn',@obj.WindowScrollWheel,'WindowButtonDownFcn',@obj.WindowButtonDown,...
                    'WindowButtonUpFcn',@obj.WindowButtonUp,'SizeChangedFcn',@obj.WindowResize,'CloseRequestFcn',@obj.WindowClose);
                obj.handle.Pointer = 'arrow';
            end
            drawnow();
        end
        %% Add an offspring
        function addOffspring(obj,h)
            h.figureNo     = length(obj.offsprings) + 1;
            obj.offsprings = [obj.offsprings,h];
            obj.offState   = [obj.offState,h.state];
            obj.offPos     = [obj.offPos;h.position];
            obj.offMoved   = [obj.offMoved,h.moved];
            obj.offPressed = [obj.offPressed,h.pressed];
        end
        %% Delete an offspring
        function delOffspring(obj,h)
            obj.offsprings(h.figureNo) = [];
            obj.offState(h.figureNo)   = [];
            obj.offPos(h.figureNo,:)   = [];
            obj.offMoved(h.figureNo)   = [];
            obj.offPressed(h.figureNo) = [];
            temp = num2cell(1:length(obj.offsprings));
            [obj.offsprings.figureNo]  = temp{:};
        end
        %% Mouse motion callback
        function WindowButtonMotion(obj,hObject,eventdata)
            mouse   = obj.handle.CurrentPoint;
            current = obj.offPos(:,1)<=mouse(1) & obj.offPos(:,2)<=mouse(2) & obj.offPos(:,1)+obj.offPos(:,3)>=mouse(1) & obj.offPos(:,2)+obj.offPos(:,4)>=mouse(2);
            current = current' & obj.offState;
            % The objects which the mouse is moving in
            for i = find(current & ~obj.offMoved)
                obj.offsprings(i).moved = true;
                obj.offsprings(i).moveIn();
            end
            % The objects which the mouse is moving out
            for i = find(~current & obj.offMoved)
                obj.offsprings(i).moved = false;
                if ~obj.offsprings(i).pressed
                    obj.offsprings(i).moveOut();
                end
            end
            % The objects which the mouse is pressing
            for i = find(obj.offPressed)
                obj.offsprings(i).mouseDrag(mouse-obj.mouseLoc);
            end
            obj.mouseLoc = mouse;
        end
        %% Mouse scrolling callback
        function WindowScrollWheel(obj,hObject,eventdata)
            for i = find(obj.offMoved)
                obj.offsprings(i).mouseScroll(eventdata.VerticalScrollCount>0);
            end
        end
        %% Mouse down callback
        function WindowButtonDown(obj,hObject,eventdata)
            if ~isempty(obj.handle.CurrentObject) && ismember(obj.handle.SelectionType,{'normal','open'})
                for i = find([obj.offsprings.handle]==obj.handle.CurrentObject & obj.offState)
                    obj.offsprings(i).pressed = true;
                    obj.offsprings(i).mouseDown();
                end
            end
        end
        %% Mouse up callback
        function WindowButtonUp(obj,hObject,eventdata)
            for i = find(obj.offPressed)
                obj.offsprings(i).pressed = false;
                if ~obj.offsprings(i).moved
                    obj.offsprings(i).moveOut();
                else
                    obj.offsprings(i).mouseUp();
                end
            end
        end
        %% Window resize callback
        function WindowResize(obj,hObject,eventdata)
            obj.position = obj.handle.Position;
        end
        %% Window close callback
        function WindowClose(obj,hObject,eventdata)
            delete(obj.handle);
            delete(obj);
        end
    end
end