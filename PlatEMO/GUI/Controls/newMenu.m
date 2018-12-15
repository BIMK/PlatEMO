classdef newMenu < newPanel
%newMenu - Menu.
%
%   h = newMenu(Parent,Pos) creates a menu which has a parent of Parent and
%   a position of Pos.
%
%   Use h.add(Str) or h.del(Index) to add an item with text Str, or delete
%   several items.
%
%   If the text of an item is empty, the item will be a line instead of an
%   option.
%
%   Use h.items(i) to obtain the i-th item in h.
%
%   If h.items(i).value = true and h.items(i).choosed = true, an icon will
%   be shown on the left of the i-th item.
%
%   The menu is invisible after being created, use h.show() to make it
%   visible.
%
%   Example:
%       newMenu(f,[10 10 100 100,1 0 1 0])
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

    properties(SetAccess = protected)
        items;       	% All the option items
        valid;          % Whether each item can be used
    end
    methods
        %% Constructor
        function obj = newMenu(parent,pos)
            obj@newPanel(parent.figure,[1 1 pos(3) 20,1 0 1 0],[],'visible',false);
        end
        %% Add one item
        function add(obj,icon,str,varargin)
            if ~isempty(str)
                obj.items = [obj.items,newMenu_Item(obj,[1 1 obj.position(3) 20,1 1 0 1],icon,str,'HorizontalAlignment','left','FontSize',10,'ForegroundColor',[.3 .3 .3],varargin{:})];
            else
                obj.items = [obj.items,newLine(obj,[30 1 obj.position(3)-30 1,1 1 0 1])];
                obj.items(end).handle.HighlightColor = [.7 .7 .7];
            end
            obj.updateLocation();
        end
        %% Delete a number of items
        function del(obj,Index)
            delete(obj.items(Index));
            obj.items(Index) = [];
            obj.updateLocation();
        end
        %% Show the menu
        function show(obj)
            % If the container of the caller object is above the menu, the
            % menu will be covered so that cannot be seen
            if obj.figure.handle.CurrentPoint(1) + (obj.position(3)-10) < obj.figure.position(3)
                Position(1) = obj.figure.handle.CurrentPoint(1) - 10;
            else
                Position(1) = obj.figure.handle.CurrentPoint(1) - (obj.position(3)-10);
            end
            if obj.figure.handle.CurrentPoint(2) - (obj.position(4)-10) > 0
                Position(2) = obj.figure.handle.CurrentPoint(2) - (obj.position(4)-10);
            else
                Position(2) = obj.figure.handle.CurrentPoint(2) - 10;
            end
            obj.position(1:2) = Position(1:2);
            obj.visible = true;
        end
    end
    methods(Access = ?newGUI)
        function moveOut(obj)
            obj.visible = false;
        end
    end
    methods(Access = protected)
        %% Update the location of each items
        function updateLocation(obj)
            if isempty(obj.items)
                obj.position(4) = 20;
            else
                obj.valid       = arrayfun(@(s)~isa(s,'newLine'),obj.items);
                obj.position(4) = sum(obj.valid)*22 + sum(~obj.valid)*5 + 4;
                start           = obj.position(4) - 1;
                for i = 1 : length(obj.items)
                    if obj.valid(i)
                        obj.items(i).position(2) = start - 22;
                        start = start - 22;
                    else
                        obj.items(i).position(2) = start - 4;
                        start = start - 5;
                    end
                end
            end
        end
    end
end