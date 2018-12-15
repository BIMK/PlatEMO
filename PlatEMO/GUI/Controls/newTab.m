classdef newTab < newPanel
%newTab - Tabgroup.
%
%   h = newTab(Parent,Pos,Str) creates a panel which has a parent of Parent
%   and a position of Pos. The panel contains multiple tab panels, where
%   Str is the set of names of the panels.
%
%   Use h.panels(i) to obtain the i-th tab panel in h.
%
%   Example:
%       newTab(f,[10 10 100 100,1 0 1 0],{'tab1','tab2','tab3'})
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
        value = 1;     	% The index of panel being chosen
    end
    properties(SetAccess = protected)
        buttons;        % All the tab buttons
        panels;         % All the tab panels
    end
    methods
        %% Constructor
        function obj = newTab(parent,pos,Str)
            obj@newPanel(parent,pos,parent.color);
            % Create the sub panels
            for i = 1 : length(Str)
                obj.panels  = [obj.panels,newPanel(obj,[1 1 pos(3) pos(4)-29,1 1 0 1],[])];
            end
            % Create the top buttons
            loc = 1;
            for i = 1 : length(Str)
                obj.buttons = [obj.buttons,newTab_Button(obj,[loc,pos(4)-29,30+length(Str{i})*10,30,1 0 0 1],Str{i},'callback',@(~,~)obj.set('value',i))];
                loc = loc + obj.buttons(i).position(3);
            end
            obj.value = 1;
        end
        %% The set method of obj.value
        function set.value(obj,value)
            % Update the current panel by obj.value
            obj.value = min(max(round(value),1),length(obj.buttons));
            [obj.buttons([1:obj.value-1,obj.value+1:end]).value]  = deal(false);
            [obj.panels([1:obj.value-1,obj.value+1:end]).visible] = deal(false);
            obj.buttons(obj.value).value  = true;
            obj.panels(obj.value).visible = true;
        end
    end
end