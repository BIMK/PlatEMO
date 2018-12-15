classdef newEdit < newGUI
%newEdit - Edit box.
%
%	h = newEdit(Parent,Pos,...) creates an edit box which has a parent of
%	Parent and a position of Pos.
%
%	Use h.string to get the text of the edit box.
%
%	Example:
%       newEdit(f,[10 10 100 100,1 0 1 0])
%
%	See also newGUI
 
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        %% Constructor
        function obj = newEdit(parent,pos,varargin)
            handle = uicontrol('FontName','Microsoft YaHei','FontSize',11,...
                               'Parent',parent.handle,'Style','edit','Units','pixels','Position',pos(1:4));
            obj@newGUI(3,parent,handle,pos(5:8),varargin{:});
            obj.handle.set('Callback',obj.callback);
            % Update obj.handle.Enable
            obj.newListener(obj,'state','PostSet',@obj.updateEnable);
        end
        %% Get the real string of the object
        function str = string(obj)
            str = obj.handle.String;
        end
    end
    methods(Access = protected)
        %% Update obj.handle.Enable
        function updateEnable(obj,hObject,eventdata)
            if obj.state
                obj.handle.Enable = 'on';
            else
                obj.handle.Enable = 'off';
            end
        end
    end
end