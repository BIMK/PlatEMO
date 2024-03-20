classdef uicontext < handle
%uicontext - Context menu with push buttons.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        fig;        % The figure
        panel;      % The main panel
        items;      % The menu items
        gaps;       % Labels for covering the borders of items
        listener;	% Listener
    end
    methods
        %% Constructor
        function obj = uicontext(fig,width)
            obj.fig   = fig;
            obj.panel = uipanel(obj.fig,'Position',[0 0 width 0],'Visible',false);
        end
        %% Add a new item
        function add(obj,str,icon,cb)
        	obj.items = [obj.items,uibutton(obj.panel,'Position',[-2 0 obj.panel.Position(3)+5 25],'Text',str,'HorizontalAlignment','left','BackgroundColor',[.94 .94 .94],'Icon',icon,'ButtonPushedFcn',cb)];
            obj.gaps  = [obj.gaps,uipanel(obj.panel,'Position',[-2 0 obj.panel.Position(3)+5 3],'BorderType','none')];
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
                if obj.fig.CurrentPoint(1)-10+obj.panel.Position(3) < obj.fig.Position(3)
                    obj.panel.Position(1) = obj.fig.CurrentPoint(1) - 10;
                else
                    obj.panel.Position(1) = obj.fig.CurrentPoint(1) + 10 - obj.panel.Position(3);
                end
                if obj.fig.CurrentPoint(2)+10-obj.panel.Position(4) > 0
                    obj.panel.Position(2) = obj.fig.CurrentPoint(2) + 10 - obj.panel.Position(4);
                else
                    obj.panel.Position(2) = obj.fig.CurrentPoint(2) - 10;
                end
                obj.panel.Visible = true;
                obj.listener      = obj.fig.addlistener('CurrentPoint','PostSet',@obj.cb_motion);
            end
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
    end
end