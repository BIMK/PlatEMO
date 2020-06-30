classdef newPopmenu < newPanel
%newPopmenu - Popup menu.
%
%   h = newPopmenu(Parent,Pos,Color,Str,...) creates a multi-column popup
%   menu which has a parent of Parent, a position of Pos and a main color
%   of Color. Str is the item group of the menu, Str(:,1) is the set of
%   group names, and Str(:,2) is the set of item names in each group.
%
%   h.index denotes the index of the current item of the menu. h.callback
%   will be invoked when h.index is set.
%
%   This kind of popup menu is non-extensible, if you want to create a
%   popup menu whose items can be added or deleted, please use newPopmenu2.
%
%   Example:
%       newPopmenu(f,[10 10 100 100,1 0 1 0],[0 0 1],...
%                  {'Group1',{'item1','item2'};'Group2',{'item1','item2'}})
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
        index = 0;                  % The index of the selected item
    end
    properties(SetAccess = ?newGUI, SetObservable)
        posChanged = false;         % Whether need to update the position of the object
    end
    properties(SetAccess = protected)
        button;                     % The button for the popup menu
        label;                      % The string of the selected item
        panel;                      % The insider panel 
        len;                        % Number of items in each group
        labels;                     % All the labels
        items;                      % All the option items
        slider;                     % The slider
    end
    methods
        %% Constructor
        function obj = newPopmenu(parent,pos,color,strGroup,varargin)
            % Create the button
            button = newPanel(parent,[pos(1:3) 20,pos(5:8)],[]);
            % Create the main panel
            obj@newPanel(parent.figure,[button.absPos(1),button.absPos(2)-200,button.absPos(3)*2,200,1 0 1 0],[],'visible',false,varargin{:});
            obj.button = button;
            % Creat the label and popup label on the button
            obj.label = newLabel(obj.button,[1 1 pos(3) 19,1 1 1 0],'','ForegroundColor',color,'callback',@obj.cb_show);
            newButtonSpecial(obj.button,[pos(3)-20 1 19 19,0 1 1 0],1,obj.color,'callback',@obj.cb_show);
            % Create the inside panel
            obj.panel = newPanel(obj,[1 1 obj.position(3:4),1 1 1 1],parent.color);
            % Create the slider
            temp = newPanel(obj,[obj.panel.position(3)-5,0,5,obj.panel.position(4),0 1 1 1],[]);
            obj.slider = newLabel(temp,[1 1 temp.position(3) temp.position(4)/3,1 1 1 0],'','BackgroundColor',color);
            % Create the labels and items
            obj.len = cellfun(@(S)length(S),strGroup(:,2));
            for i = 1 : size(strGroup,1)
                obj.labels = [obj.labels,newLabel(obj.panel,[1 20*i obj.panel.position(3) 15,1 1 0 1],[' ',strGroup{i,1}],...
                                                  'HorizontalAlignment','left','BackgroundColor',color,'FontSize',9,'ForegroundColor',[1 1 1])];
                for j = 1 : length(strGroup{i,2})
                    k = length(obj.items);
                    if length(strGroup{i,2}{j}) < 9
                        FontSize = 11;
                    else
                        FontSize = 10;
                    end
                    obj.items = [obj.items,newLabelButton(obj.panel,[1 1 115 20,1 0 1 0],strGroup{i,2}{j},'ForegroundColor',color,'FontSize',FontSize,'callback',@(~,~)obj.set('index',k+1))];
                end
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
                value = '';
            end
        end
    end
    methods(Access = ?newGUI)
        function moveOut(obj)
            obj.visible = false;
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
        %% The callback of obj.label
        function cb_show(obj,hObject,eventdata)
            % Update the position of obj
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
            obj.position([1,3]) = [obj.button.absPos(1),obj.button.absPos(3)*2];
            nColumn             = max(floor(obj.panel.position(3)/115),1);
            high                = 20*length(obj.len) + sum(ceil(obj.len/nColumn)*25);
            obj.position([2,4]) = [obj.button.absPos(2)-min(high,200),min(high,200)];
            % Calculate the new position of obj.panel
            obj.panel.position([2,4]) = [obj.position(4)-max(obj.position(4),high)+1,max(obj.position(4),high)];
            % Calculate the new position of labels and items
            start  = obj.panel.position(4);
            cumlen = [0;cumsum(obj.len)];
            for i = 1 : length(obj.len)
                obj.labels(i).position(2) = start - 15;
                for j = 1 : obj.len(i)
                    obj.items(cumlen(i)+j).position(1:2) = [mod(j-1,nColumn)*115+1,start-15-ceil(j/nColumn)*25];
                end
                start = obj.items(cumlen(i+1)).position(2) - 5;
            end
            % Calculate the new position of slider
            obj.slider.position(2)    = obj.position(4) - obj.slider.position(4);
            obj.slider.parent.visible = obj.position(4) < obj.panel.position(4);
        end
    end
end