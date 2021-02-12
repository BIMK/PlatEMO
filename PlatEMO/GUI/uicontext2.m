classdef uicontext2 < handle
%uicontext2 - Context menu with state buttons.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        value;      % Index of the pressed button
        string;     % String of selected item
    end
    properties(SetAccess = private)
        fig;        % The figure
        panel;      % The main panel
        items;      % The menu items
        gaps;       % Labels for covering the borders of items
        callback;   % Value changed callback function
        listener;	% Listener
    end
    methods
        %% Constructor
        function obj = uicontext2(fig,cb)
            obj.fig      = fig;
            obj.panel    = uipanel(obj.fig,'Position',[0 0 110 0],'Visible',false);
            obj.callback = cb;
        end
        %% Add a new item
        function add(obj,str,gap)
        	obj.items = [obj.items,uibutton(obj.panel,'Position',[-3 0 120 25],'Text',['    ',str],'HorizontalAlignment','left','BackgroundColor',[.94 .94 .94],'ButtonPushedFcn',{@obj.cb_button,length(obj.items)+1})];
            obj.gaps  = [obj.gaps,uipanel(obj.panel,'Position',[-3 0 120 (1-gap)*2+1],'BorderType','none')];
            obj.panel.Children = obj.panel.Children([1,3:end,2]);
        end
        %% Flush the menu
        function flush(obj)
            if ~isempty(obj.items)
                obj.panel.Position(4) = 23*length(obj.items) - 2;
                loc = obj.panel.Position(4) + 1;
                for i = 1 : length(obj.items)
                    loc = loc - 23;
                    obj.items(i).Position(2) = loc;
                    obj.gaps(i).Position(2)  = loc - 1;
                end
            end
        end
        %% Show the menu
        function show(obj)
            if ~isempty(obj.items)
                obj.panel.Position(1:2) = [obj.fig.CurrentPoint(1)-10,obj.fig.CurrentPoint(2)+10-obj.panel.Position(4)];
                obj.panel.Visible       = true;
                obj.listener            = obj.fig.addlistener('CurrentPoint','PostSet',@obj.cb_motion);
            end
        end
        %% Change obj.value
        function set.value(obj,p)
            [obj.items.BackgroundColor]  = deal([.94 .94 .94]);
            obj.items(p).BackgroundColor = [.6 .75 .9];
            obj.value = p;
        end
        %% Get the string of selected item
        function str = get.string(obj)
            str = strtrim(obj.items(obj.value).Text);
        end
    end
    methods(Access = private)
        %% Hide the menu when moving out of it
        function cb_motion(obj,~,~)
            if obj.fig.CurrentPoint(1)<obj.panel.Position(1) || obj.fig.CurrentPoint(1)>sum(obj.panel.Position([1,3])) || obj.fig.CurrentPoint(2)<obj.panel.Position(2) || obj.fig.CurrentPoint(2)>sum(obj.panel.Position([2,4]))
                obj.panel.Visible = false;
                delete(obj.listener);
            end
        end
        %% Change the selected button and execute the callback
        function cb_button(obj,~,~,p)
            obj.value = p;
            obj.panel.Visible = false;
        	delete(obj.listener);
            if isa(obj.callback,'function_handle')
                obj.callback([],[]);
            else
                obj.callback{1}([],[],obj.callback{2:end});
            end
        end
    end
end