classdef newGUI < handle & matlab.mixin.Heterogeneous
%newGUI - The superclass of all the GUI classes.
%
%   This is the superclass of all the other GUI classes. This class cannot
%   be instantiated.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetObservable)
        enable  = true;                 % Whether the object is enabled
        visible = true;                 % Whether the object is visible
        state   = true;                 % Whether the object is active
        moved   = false;                % Whether the object is being moved over
        pressed = false;                % Whether the object is being pressed
        position;                       % Position of the object
        callback     = @(varargin)[];   % Function invoked when click the control
        movecallback = @(varargin)[];   % Function invoked when move in the control
    end
    properties(SetAccess = ?newGUI, SetObservable)
        absPos;                         % Absolute position of the object (transitive)
        oldAbsPos;                      % Last value of obj.absPos
        oldParentAbsPos;                % Last absolute position of its parent when it is visible
        children;                       % All the children of this object
        figureNo;                       % The No. of this object in the figure
    end
    properties(SetAccess = private)
        objtype;                        % Type of object(1.figure 2.container 3.control)
        parent;                         % The parent object
        handle;                         % Main GUI object
        figure;                         % The figure object
        sizeLock;                       % The size of the object locked to its parent (left,right,bottom,top)
        absOffset = -1;                 % Offset for absPos calculation
        listeners;                      % All the listeners in the object
    end
    methods(Access = protected)
        %% Constructor
        function obj = newGUI(objtype,parent,handle,sizeLock,varargin)
            % Initialise the basic properties
            obj.objtype  = objtype;
            obj.parent   = parent;
            obj.handle   = handle;
            obj.sizeLock = [sizeLock(1)*10+sizeLock(2),sizeLock(3)*10+sizeLock(4)];
            % Add the object to its parent and figure
            if obj.objtype == 1
                obj.figure = obj;
            else
                obj.figure = obj.parent.figure;
                obj.figure.addOffspring(obj);
                obj.parent.children = [obj.parent.children,obj];
            end
            % Initialise the properties
            if ~isempty(varargin)
                IsString = find(cellfun(@ischar,varargin(1:end-1)));
                Loc      = IsString(ismember(varargin(IsString),properties(obj)));
                obj.set(varargin{[Loc;Loc+1]});
                varargin([Loc;Loc+1]) = [];
                if ~isempty(varargin)
                    set(obj.handle,varargin{:});
                end
            end
            % Set other properties
            obj.position = obj.handle.Position;
            obj.state    = obj.enable & obj.visible & obj.parent.state;
            if obj.objtype == 1
                obj.absPos = [0,0,obj.position(3:4)];
                obj.oldParentAbsPos = 0;
            else
                obj.absPos = obj.position + [obj.parent.absPos(1:2)+obj.absOffset,0,0];
                obj.oldParentAbsPos = obj.parent.oldAbsPos;
            end
            obj.oldAbsPos = obj.absPos;
        end
        %% Set the property
        function set(obj,varargin)
            for i = 1 : 2 : length(varargin)
                obj.(varargin{i}) = varargin{i+1};
            end
        end
    end
    methods
        %% Destructor
        function delete(obj)
            try
                % Delete the object from the figure
                if isvalid(obj.figure)
                    obj.figure.delOffspring(obj);
                end
                % Delete the object from its parent
                if isa(obj.parent,'newGUI') && isvalid(obj.parent)
                    obj.parent.children([obj.parent.children.handle]==obj.handle) = [];
                end
                % Delete the handle objects and the newGUI objects of
                % its children
                Values = cellfun(@(S)obj.(S),properties(obj),'UniformOutput',false);
                Del    = cellfun(@(S)isa(S,'handle'),Values) & ~ismember(properties(obj),{'parent','figure'});
                cellfun(@(S)delete(S),Values(Del));
            catch
            end
        end
        %% Add a new listener
        function newListener(obj,varargin)
            obj.listeners = [obj.listeners,addlistener(varargin{:})];
        end
        %% Get the background color of the object
        function value = color(obj)
            value = obj.handle.BackgroundColor;
        end
        %% Move the object to the top of its parent's children
        function movetop(obj)
            index = find([obj.parent.handle.Children]==obj.handle,1);
            obj.parent.handle.Children = obj.parent.handle.Children([index,1:index-1,index+1:end]);
        end
        %% All the set and get methods
        function set.enable(obj,value)
            obj.enable = value;
            % Update obj.state
            obj.state = obj.enable & obj.visible & obj.parent.state;
        end
        function set.visible(obj,value)
            oldVisible  = obj.visible;
            obj.visible = value;
            % Update obj.state
            obj.state = obj.enable & obj.visible & obj.parent.state;
            if value && ~oldVisible
                % Update obj.handle.Visible
                obj.handle.Visible = 'on';
                % Update obj.position
                Position = obj.position;
                switch obj.sizeLock(1)
                    case 1
                        Position(1) = Position(1) + obj.parent.oldAbsPos(3) - obj.oldParentAbsPos(3);
                    case 0
                        Position([1,3]) = Position([1,3]).*obj.parent.oldAbsPos(3)./obj.oldParentAbsPos(3);
                    case 11
                        Position(3) = Position(3) + obj.parent.oldAbsPos(3) - obj.oldParentAbsPos(3);
                end
                switch obj.sizeLock(2)
                    case 1
                        Position(2) = Position(2) + obj.parent.oldAbsPos(4) - obj.oldParentAbsPos(4);
                    case 0
                        Position([2,4]) = Position([2,4]).*obj.parent.oldAbsPos(4)./obj.oldParentAbsPos(4);
                    case 11
                        Position(4) = Position(4) + obj.parent.oldAbsPos(4) - obj.oldParentAbsPos(4);
                end
                obj.position = Position;
            elseif ~value && oldVisible
                % Update obj.handle.Visible
                obj.handle.Visible = 'off';
                % Update obj.oldParentAbsPos
            	obj.oldParentAbsPos = obj.parent.oldAbsPos;
            end
        end
        function set.state(obj,value)
            obj.state = value;
            % Update its state in the figure
            if obj.objtype > 1
                obj.figure.offState(obj.figureNo) = value;
            end
            % Update obj.children.state
            if obj.objtype < 3
                obj.updateChildrenState();
            end
        end
        function set.moved(obj,value)
            obj.moved = value;
            obj.figure.offMoved(obj.figureNo) = value;
        end
        function set.pressed(obj,value)
            obj.pressed = value;
            obj.figure.offPressed(obj.figureNo) = value;
        end
        function set.position(obj,value)
            obj.position = value;
            % Update obj.handle.Position
            obj.handle.Position = [value(1:2),max(1,value(3:4))];
            % Update obj.absPos
            if obj.objtype == 1
                obj.absPos = [0,0,value(3:4)];
            else
                obj.absPos = value + [obj.parent.absPos(1:2)+obj.absOffset,0,0];
                obj.figure.offPos(obj.figureNo,:) = obj.absPos;
            end
            % Update obj.children.position
            if obj.objtype < 3
                obj.updateChildrenPosition();
            end
        end
    end
    methods(Access = ?newGUI)
        %% Update children's state
        function updateChildrenState(obj)
            if ~isempty(obj.children)
                State = num2cell([obj.children.enable]&[obj.children.visible]&obj.state);
                [obj.children.state] = State{:};
            end
        end
        %% Update children's position
        function updateChildrenPosition(obj)
            if ~isempty(obj.children) && any(obj.oldAbsPos~=obj.absPos)
                Children = obj.children([obj.children.visible]);
                if ~isempty(Children)
                    if obj.oldAbsPos(3:4) == obj.absPos(3:4)
                        % Only update obj.children.absPos
                        AbsPos = cat(1,Children.absPos)+repmat(obj.absPos-obj.oldAbsPos,length(Children),1);
                        obj.figure.offPos([Children.figureNo],:) = AbsPos;
                        AbsPos = num2cell(AbsPos,2);
                        [Children.absPos] = AbsPos{:};
                        arrayfun(@(h)h.updateChildrenPosition(),Children([Children.objtype]<3));
                    else
                        % Update obj.children.position
                        for h = Children
                            Position = h.position;
                            switch h.sizeLock(1)
                                case 1
                                    Position(1) = Position(1) + obj.absPos(3) - obj.oldAbsPos(3);
                                case 0
                                    Position([1,3]) = Position([1,3]).*obj.absPos(3)./obj.oldAbsPos(3);
                                case 11
                                    Position(3) = Position(3) + obj.absPos(3) - obj.oldAbsPos(3);
                            end
                            switch h.sizeLock(2)
                                case 1
                                    Position(2) = Position(2) + obj.absPos(4) - obj.oldAbsPos(4);
                                case 0
                                    Position([2,4]) = Position([2,4]).*obj.absPos(4)./obj.oldAbsPos(4);
                                case 11
                                    Position(4) = Position(4) + obj.absPos(4) - obj.oldAbsPos(4);
                            end
                            h.position = Position;
                        end
                    end
                end
            end
            % Update obj.oldAbsPos
            obj.oldAbsPos = obj.absPos;
        end
        %% All the inner callback functions
        function moveIn(obj)
            if obj.objtype >= 3
                obj.movecallback(obj,[]);
            end
        end
        function moveOut(obj)
        end
        function mouseScroll(obj,moveDown)
        end
        function mouseDown(obj)
        end
        function mouseUp(obj)
            if obj.objtype >= 3
                obj.callback(obj,[]);
            end
        end
        function mouseDrag(obj,offset)
        end
    end
end