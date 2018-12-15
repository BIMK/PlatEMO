classdef newIcon < newGUI
%newIcon - Static icon.
%
%   h = newIcon(Parent,Pos,Icon,...) creates a static icon which has a
%   parent of Parent, a position of Pos and an icon of Icon.
%
%   Example:
%       newIcon(f,[10 10 100 100,1 0 1 0],cdata)
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

    properties(SetAccess = protected)
        icon;               % The icon data
        CDataNormal;        % The CData of normal state
        CDataUnable;        % The CData of unable state
    end
    methods
        %% Constructor
        function obj = newIcon(parent,pos,icon)
            % Create the button-like icon
            handle = uicontrol('Parent',parent.handle,'Style','pushbutton','Enable','inactive','Units','pixels','Position',pos(1:4));
            obj@newGUI(3,parent,handle,pos(5:8));
            obj.icon = icon;
            % Update the button style by obj.moved, obj.pressed, obj.value
            % or obj.state
            obj.newListener(obj,'state','PostSet',@obj.updateStyle);
            % Update the CData by obj.position
            obj.newListener(obj.handle,'SizeChanged',@obj.updateCData);
            obj.updateCData();
        end
    end
    methods(Access = protected)
    	%% Update the CData
        function updateCData(obj,hObject,eventdata)
            try
                pos = floor(obj.handle.Position)+2;
                % Expand the icon data
                CData = nan(pos(4),pos(3),3);
                sizeI = floor([size(obj.icon,1),size(obj.icon,2)]/2);
                CData(ceil(pos(4)*0.50)+(-sizeI(1)+1:size(obj.icon,1)-sizeI(1)),...
                      ceil(pos(3)*0.50)+(-sizeI(2):size(obj.icon,2)-sizeI(2)-1),:) = obj.icon;
                transparent = repmat(any(isnan(CData),3),1,1,3);
                len         = sum(sum(sum(transparent)))/3;
                CData       = min(max(CData,0),1);
                % The foreground color
                parentColor     = repmat(reshape(obj.parent.color,1,1,3),pos(4),pos(3),1);
                obj.CDataNormal = CData;
                obj.CDataUnable = (CData+parentColor)/2;
                % The background color
                obj.CDataNormal(transparent) = repmat(obj.parent.color,len,1);
                obj.CDataUnable(transparent) = repmat(obj.parent.color,len,1);
                obj.updateStyle();
            catch
            end
        end
        %% Update the button style
        function updateStyle(obj,hObject,eventdata)
            if obj.state
            	obj.handle.CData = obj.CDataNormal;
            else
            	obj.handle.CData = obj.CDataUnable;
            end
        end
    end
end