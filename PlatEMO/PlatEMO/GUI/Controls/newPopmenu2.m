classdef newPopmenu2 < newPanel
%newPopmenu2 - Popup menu II.
%
%   h = newPopmenu2(Parent,Pos,Color,Del,Str,...) creates a popup menu
%   which has a parent of Parent, a position of Pos and a main color of
%   Color. Del indicates whether there is a delete button. Str is the set
%   of item names.
%
%   h.index denotes the index of the current item of the menu. h.callback
%   will be invoked when h.index is set. h.delcallback will be invoked when
%   delete items from the menu.
%
%   Use h.add(Str) or h.del(Index) to add items with text Str, or delete
%   several items.
%
%   Example:
%       newPopmenu2(f,[10 10 100 100,1 0 1 0],[0 0 1],1,...
%                   {'item1','item2','item3'})
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
        index = 0;                      % The index of the selected item
        delcallback = @(varargin)[];	% Callback of deleting items
    end
	properties(SetAccess = ?newGUI, SetObservable)
        posChanged = false;             % Whether need to update the position of the object
    end
    properties(SetAccess = protected)
        button;                         % The button for the popup menu
        label;                          % The string of the selected item
        panel;                          % The selection panel (insider)
        items;                          % All the option items
        delButton;                      % The delete button
        slider;                         % The slider
    end
    methods
        %% Constructor
        function obj = newPopmenu2(parent,pos,color,del,Str,varargin)
            % Create the button
            button = newPanel(parent,[pos(1:3) 20,pos(5:8)],[]);
            % Create the main panel
            obj@newPanel(parent.figure,[button.absPos(1),button.absPos(2)-100,button.absPos(3),100,1 0 1 0],[],'visible',false,varargin{:});
            obj.button = button;
            % Create the label and popup label on the button
            obj.label = newLabel(obj.button,[1 1 pos(3) 19,1 1 1 0],'','ForegroundColor',color,'callback',@obj.cb_show);
            newButtonSpecial(obj.button,[pos(3)-20 1 19 19,0 1 1 0],1,obj.color,'callback',@obj.cb_show);
            % Create the inside panel
            obj.panel = newPanel(obj,[1 1 obj.position(3:4),1 1 1 1],parent.color);
            % Create the slider
            temp = newPanel(obj,[obj.panel.position(3)-5,0,5,obj.panel.position(4),0 1 1 1],[]);
            obj.slider = newLabel(temp,[1 1 temp.position(3) temp.position(4)/3,1 1 1 0],'','BackgroundColor',color);
            % Create the items
            if nargin > 4 && ~isempty(Str)
                obj.add(Str);
            end
            % Create the delete button
            if del
                obj.delButton = newButtonSpecial(obj.panel,[obj.panel.position(3)-22 1 20 20,0 1 1 0],3,(obj.color+[.8 .8 1])/2,'ForegroundColor',[.6 .6 1]);
                obj.delButton.visible = false;
            end
            % Change the position of obj according to obj.button.position
            obj.newListener(obj.button,'position','PostSet',@(~,~)obj.set('posChanged',true));
            obj.newListener(obj.button,'moved','PostSet',@(~,~)obj.set('visible',obj.visible&(obj.moved|obj.button.moved)));
            obj.updateLocation();
        end
        %% The set method of obj.index
        function set.index(obj,value)
            obj.index = value;
            % Hide the object
            obj.visible = false;
            % Update the string of popup menu
            if obj.index
                obj.label.handle.String = obj.string;
                set([obj.items.handle],'FontWeight','normal');
                set(obj.items(obj.index).handle,'FontWeight','bold');
            else
                obj.label.handle.String = '';
            end
            obj.callback(obj,obj.index);
        end
        %% Get the string of the selected item
        function value = string(obj)
            if obj.index > 0
                value = obj.items(obj.index).handle.String;
            else
                value = 0;
            end
        end
        %% Add items
        function add(obj,Str)
            if ~iscell(Str)
                Str = {Str};
            end
            % Add an item
            for i = 1 : length(Str)
                obj.items = [obj.items,newPopmenu2_Item(obj.panel,[1 1 obj.panel.position(3) 20,1 1 1 0],Str{i},'ForegroundColor',obj.label.foregroundcolor)];
            end
            % Update the location of objects
            obj.updateLocation();
             % Move the delete button to the top
            if ~isempty(obj.delButton)
                obj.delButton.movetop();
            end
        end
        %% Delete items
        function del(obj,Index)
            obj.delcallback(obj,Index);
            % Delete the items
            delete(obj.items(Index));
            obj.items(Index) = [];
            % Update obj.index
            obj.index = min(obj.index-sum(Index<obj.index),length(obj.items));
            % Update the location of objects
            obj.updateLocation();
        end
    end
    methods(Access = ?newGUI)
        function moveOut(obj)
            obj.visible = false;
            if ~isempty(obj.delButton)
                obj.delButton.visible = false;
            end
        end
        function mouseScroll(obj,moveDown)
            if obj.panel.position(4) > obj.position(4)
                % Update the vertical location of obj.panel
                if moveDown
                    obj.panel.position(2) = max(min(obj.panel.position(2)+30,1),obj.position(4)-obj.panel.position(4)+1);
                else
                    obj.panel.position(2) = max(min(obj.panel.position(2)-30,1),obj.position(4)-obj.panel.position(4)+1);
                end
                % Update the vertical location of obj.slider
                obj.slider.position(2) = -(1-obj.panel.position(2)) / (-obj.position(4)+obj.panel.position(4)) * (1-obj.position(4)+obj.slider.position(4));
            end
        end
    end
    methods(Access = protected)
        %% The callback of obj.label and obj.button
        function cb_show(obj,hOject,eventdata)
            % Update the position of obj.list
            if obj.posChanged
                obj.posChanged = false;
                obj.updateLocation();
            end
            % Show the list or not
            obj.visible = ~obj.visible;
        end
        %% Update the location of obj.handle, panel, labels and items by the position of popmenu
        function updateLocation(obj)
            % Calculate the new position of obj.handle
            high         = length(obj.items)*22 + 4;
            obj.position = [obj.button.absPos(1),obj.button.absPos(2)-max(min(high,100),10),obj.button.absPos(3),max(min(high,100),10)];
            % Calculate the new position of obj.panel
            obj.panel.position([2,4]) = [obj.position(4)-max(obj.position(4),high)+1,max(obj.position(4),high)];
            % Calculate the new position of labels and items
            for i = 1 : length(obj.items)
                obj.items(i).position(2) = obj.panel.position(4) - 1 - i*22;
            end
            % Calculate the new position of slider
            obj.slider.position(2)    = obj.position(4) - obj.slider.position(4);
            obj.slider.parent.visible = obj.position(4) < obj.panel.position(4);
        end
    end
end