classdef ParameterList < newPanel
%ParameterList - Panel for parameter setting.
%
%   h = ParameterList(Parent,Pos,Note,Icon...) creates a panel for
%   parameter setting which has a parent of Parent and a position of Pos.
%   Note is the corresponding note label of h, and Icon is the set of icons
%   for the menu.
%
%   Use h.update(Operation,Str,Color,Type,Replace) to add or replace the
%   items in h.
%
%   Example:
%       ParameterList(f,[10 10 100 100,1 0 1 0],Obj)
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

    properties(SetAccess = private)
        slider;                     % The main slider
        items;                      % The setting groups
        menus;                      % The menu objects
        note;                       % The corresponding note label object
    end
    properties(SetAccess = ?ParameterList_Item)
        currentIndex;               % The index of item which calls the menu
    end
    methods
        function obj = ParameterList(parent,pos,note,Icon,varargin)
            obj@newPanel(parent,pos,parent.color,varargin{:});
            obj.note = note;
            % Create the slider
            obj.slider = newSlider(obj,[pos(3)-17,1,18,pos(4),0 1 1 1],'visible',false,'callback',@obj.cb_updateList);
            % Create the menus
            menu1 = newMenu(parent,[1 1 120 1,1 0 1 0]);
            menu1.add(Icon.openfile,'Open file','callback',@obj.cb_openfile);
            menu1.add(Icon.openfolder,'Open folder','callback',@obj.cb_openfolder);
            menu1.add(Icon.scholar,'Search online','callback',@obj.cb_search);
            menu2 = newMenu(parent,[1 1 120 1,1 0 1 0]);
            menu2.add(Icon.moveup,'Move up','callback',@obj.cb_moveup);
            menu2.add(Icon.movedown,'Move down','callback',@obj.cb_movedown);
            menu2.add(Icon.delete,'Delete','callback',@obj.cb_delete);
            menu2.add([],'');
            menu2.add(Icon.openfile,'Open file','callback',@obj.cb_openfile);
            menu2.add(Icon.openfolder,'Open folder','callback',@obj.cb_openfolder);
            menu2.add(Icon.scholar,'Search online','callback',@obj.cb_search);
            obj.menus = [menu1,menu2];
            % Show or hide the slider according to obj.position
            obj.newListener(obj,'position','PostSet',@obj.cb_updateList);
        end
        %% Add or replace items
        function update(obj,operation,itemNames,color,type,replace)
            if ~iscell(itemNames)
                itemNames = {itemNames};
            else
                itemNames = reshape(unique(itemNames,'stable'),1,[]);
            end
            % Identify the new items in itemNames
            if isempty(obj.items)
                new = true(1,length(itemNames));
            else
                new = ~ismember(itemNames,{obj.items.name});
            end
            % Add or replace items
            switch operation
                case 'add'
                    % Add new items
                    for i = find(new)
                        obj.items = [obj.items,ParameterList_Item(obj,itemNames{i},color,type)];
                    end
                    changed = any(new);
                case 'replace'
                    % Delete old items
                    delete(obj.items(replace(replace<=length(obj.items)&new)));
                    % Replace with new items
                    for i = find(replace<=length(obj.items)&new)
                        obj.items(replace(i)) = ParameterList_Item(obj,itemNames{i},color,type);
                    end
                    % Add new items
                    for i = find(replace>length(obj.items)&new)
                        obj.items = [obj.items,ParameterList_Item(obj,itemNames{i},color,type)];
                    end
                    changed = any(new);
                case 'replaceall'
                    % Delete old items
                    if ~isempty(obj.items)
                        del = ~ismember({obj.items.name},itemNames);
                        for i = find(del)
                            delete(obj.items(i));
                        end
                        obj.items(del) = [];
                    else
                        del = [];
                    end
                    % Add new items
                    for i = find(new)
                        obj.items = [obj.items,ParameterList_Item(obj,itemNames{i},color,type)];
                    end
                    % Sort the items
                    if ~isempty(obj.items)
                        [~,rank]  = ismember(itemNames,{obj.items.name});
                        obj.items = obj.items(rank);
                        changed   = any(new) || any(del) || ~issorted(rank);
                    else
                        changed = any(new) || any(del);
                    end
            end
            if changed
                obj.cb_updateList();
                obj.callback(obj,[]);
            end
        end
        %% Update the location of each item
        function cb_updateList(obj,hObject,eventdata)
            high = 0;
            for h = obj.items
                high = high + h.position(4);
            end
            if high <= obj.position(4)
                obj.slider.visible = false;
                start = obj.position(4);
                width = obj.position(3);
            else
                obj.slider.visible = true;
                start = (high-obj.position(4))*obj.slider.value + obj.position(4);
                width = obj.position(3) - 20;
            end
            for h = obj.items
                start = start - h.position(4);
                h.position([2,3]) = [start,width];
            end
        end
    end
    methods(Access = ?newGUI)
        function mouseScroll(obj,moveDown)
            if moveDown
                obj.slider.value = obj.slider.value + 0.1;
            else
                obj.slider.value = obj.slider.value - 0.1;
            end
        end
    end
    methods(Access = protected)
        %% The callback of menu (move the item up)
        function cb_moveup(obj,hObject,eventdata)
            if obj.currentIndex > 1
                obj.items = obj.items([1:obj.currentIndex-2,obj.currentIndex,obj.currentIndex-1,obj.currentIndex+1:end]);
                obj.cb_updateList();
            end
        end
        %% The callback of menu (move the item down)
        function cb_movedown(obj,hObject,eventdata)
            if obj.currentIndex < length(obj.items)
                obj.items = obj.items([1:obj.currentIndex-1,obj.currentIndex+1,obj.currentIndex,obj.currentIndex+2:end]);
                obj.cb_updateList();
            end
        end
        %% The callback of menu (delete the item)
        function cb_delete(obj,hObject,eventdata)
            delete(obj.items(obj.currentIndex));
            obj.items(obj.currentIndex) = [];
            obj.cb_updateList();
            obj.callback(obj,[]);
        end
        %% The callback of menu (open file)
        function cb_openfile(obj,hObject,eventdata)
            web(['file://',which(obj.items(obj.currentIndex).name)],'-browser');
        end
        %% The callback of menu (open folder)
        function cb_openfolder(obj,hObject,eventdata)
            web(['file://',fileparts(which(obj.items(obj.currentIndex).name))],'-browser');
        end
        %% The callback of menu (search online)
        function cb_search(obj,hObject,eventdata)
            web(['https://scholar.google.com/scholar?q=%',strjoin(cellstr(dec2hex(double(obj.items(obj.currentIndex).ref))),'%')],'-browser');
        end
    end
end