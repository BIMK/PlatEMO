classdef newPopmenu3 < newLabelButton
%newPopmenu3 - Popup menu III.
%
%   h = newPopmenu3(Parent,Pos,Str,...) creates a button and menu based
%   popup menu which has a parent of Parent, a position of Pos. Str is the
%   set of item names.
%
%   h.index denotes the index of the current item of the menu. h.callback
%   will be invoked when h.index is set.
%
%   This kind of popup menu is non-extensible, if you want to create a
%   popup menu whose items can be added or deleted, please use newPopmenu2.
%
%   Example:
%       newPopmenu3(f,[10 10 100 100,1 0 1 0],{'item1','item2','item3'})
%
%   See also newGUI, newLabelButton

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetObservable)
        index = 1;      % The index of selected item in the menu
    end
    properties(SetAccess = protected)
        menu;           % Popup menu
        tip;            % Tip button
    end
    methods
        %% Constructor
        function obj = newPopmenu3(parent,pos,Icon,Str,varargin)
            obj@newLabelButton(parent,pos,'','FontSize',13,varargin{:});
            % Creat the popup menu
            obj.menu = newMenu(parent,[1 1 115 1,1 0 1 0]);
            % Create the popup label
            obj.tip = newTip(parent,[1 1 12 1],obj,'callback',@(~,~)obj.menu.show());
            % Creat the items in the menu
            valid = cellfun(@(S)~isempty(S),Str);
            for i = 1 : length(Str)
                obj.menu.add(Icon{i},Str{i},'choosed',true,'callback',@(~,~)obj.set('index',sum(valid(1:i))));
            end
            % Select the first item
            first = find(obj.menu.valid,1);
            obj.menu.items(first).value = true;
            obj.handle.String = regexp(obj.menu.items(first).handle.String,'(?<=&nbsp\s+).*','match','once');
        end
        %% The set method of obj.index
        function set.index(obj,value)
            obj.index = value;
            valid     = find(obj.menu.valid);
            [obj.menu.items(obj.menu.valid).value] = deal(false);
            obj.menu.items(valid(obj.index)).value = true;
            obj.handle.String = regexp(obj.menu.items(valid(obj.index)).handle.String,'(?<=&nbsp\s+).*','match','once');
            obj.callback(obj,obj.index);
        end
    end
    methods(Access = ?newGUI)
        function mouseUp(obj)
            obj.value = ~obj.value;
            obj.menu.show();
        end
    end
end