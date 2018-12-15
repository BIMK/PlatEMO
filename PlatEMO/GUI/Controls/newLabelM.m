classdef newLabelM < newGUI
%newLabelM - Multi-line label.
%
%   h = newLabelM(Parent,Pos,Str,...) creates a multi-line label which has
%   a parent of Parent, a position of Pos and a text of Str. The text can
%   have multiple lines (set as string cell).
%
%   Example:
%       newLabelM(f,[10 10 100 100,1 0 1 0])
%
%   See also newGUI

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
        function obj = newLabelM(parent,pos,str,varargin)
            handle = uicontrol('BackgroundColor',[0 0 0]+0.941,'FontSize',9,'String',str,'HorizontalAlignment','left','Max',2,...
                              'Parent',parent.handle,'Style','edit','Enable','inactive','Units','pixels','Position',pos(1:4));              
            obj@newGUI(3,parent,handle,pos(5:8),varargin{:});
            % Update the label style by obj.state
            obj.newListener(obj,'state','PostSet',@obj.updateStyle);
        end
    end
    methods(Access = protected)
        %% Update the label style
        function updateStyle(obj,hObject,eventdata)
            if obj.state
                obj.handle.ForegroundColor = [0 0 0];
            else
                obj.handle.ForegroundColor = [.5 .5 .5];
            end
        end
    end
end