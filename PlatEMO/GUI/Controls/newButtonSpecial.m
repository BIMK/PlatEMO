classdef newButtonSpecial < newButton
%newButtonSpecial - Button with special character.
%
%   h = newButtonSpecial(Parent,Pos,Type,Color, ...) creates a Button which
%   has a parent of Parent, a position of Pos, a character specified by
%   Type and a background color of Color.
%
%   If Color = [], this control looks like a button; otherwise it looks
%   like a label.
%
%   Example:
%       newButtonSpecial(f,[10 10 100 100,1 0 1 0],1)
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
        icon;        	% The icon data
        color;      	% The background color of the label
    end
    methods
        %% Constructor
        function obj = newButtonSpecial(parent,pos,type,color,varargin)
            % Create the label (button based label)
            obj@newButton(parent,pos,'',varargin{:});
            % Get the icon and the background color
            obj.icon  = obj.getIcon(type);
            obj.color = color;
            obj.updateCData();
        end
    end
    methods(Access = protected)
        %% Update the CData
        function updateCData(obj,hObject,eventdata)
            try
                pos = floor(obj.handle.Position)+2;
                % Expand the icon data
                CData = false(pos(4),pos(3));
                sizeI = floor([size(obj.icon,1),size(obj.icon,2)]/2);
                CData(ceil(pos(4)*0.50)+(-sizeI(1)+1:size(obj.icon,1)-sizeI(1)),...
                      ceil(pos(3)*0.50)+(-sizeI(2):size(obj.icon,2)-sizeI(2)-1)) = obj.icon;
                transparent = repmat(~CData,1,1,3);
                len         = sum(sum(sum(transparent)))/3;
                % Update the CData of the button
                if ~isempty(obj.color)
                    % The foreground color
                    obj.CDataNormal = repmat(reshape(obj.handle.ForegroundColor,1,1,3),pos(4),pos(3),1);
                    obj.CDataUnable = repmat(reshape((obj.handle.ForegroundColor+obj.parent.color)/2,1,1,3),pos(4),pos(3),1);
                    % The background color
                    obj.CDataNormal(transparent) = repmat(obj.color,len,1);
                    obj.CDataUnable(transparent) = repmat((obj.color+obj.parent.color)/2,len,1);
                else
                    % The foreground color
                    obj.CDataNormal  = repmat(0.2,pos(4),pos(3),3);
                    obj.CDataUnable  = repmat(0.5,pos(4),pos(3),3);
                    obj.CDataMove    = obj.CDataNormal;
                    obj.CDataPress   = obj.CDataNormal;
                    obj.CDataUnableP = obj.CDataUnable;
                    % The background color
                    obj.CDataNormal(transparent)  = repmat(obj.parent.color,len,1);
                    obj.CDataUnable(transparent)  = repmat(obj.parent.color,len,1);
                    obj.CDataMove(transparent)    = repmat((obj.parent.color+[.8 .8 1])/2,len,1);
                    obj.CDataPress(transparent)   = repmat((obj.parent.color+[.6 .6 1])/2,len,1);
                    obj.CDataUnableP(transparent) = repmat([.9 .9 .9],len,1);
                end
                obj.updateStyle();
            catch
            end
        end
        %% Update the label style
        function updateStyle(obj,hObject,eventdata)
            if obj.state
                if ~isempty(obj.color)
                    obj.handle.CData = obj.CDataNormal;
                elseif obj.value && obj.choosed || obj.pressed
                    obj.handle.CData = obj.CDataPress;
                elseif obj.moved
                    obj.handle.CData = obj.CDataMove;
                else
                	obj.handle.CData = obj.CDataNormal;
                end
            else
                if isempty(obj.color) && obj.value && obj.choosed
                    obj.handle.CData = obj.CDataUnableP;
                else
                    obj.handle.CData = obj.CDataUnable;
                end
            end
        end
        %% Get the icon data of specified character
        function icon = getIcon(obj,type)
            switch type
                case 1
                    icon = [1 1 1 1 1 1 1 1 1 1 1
                            0 1 1 1 1 1 1 1 1 1 0
                            0 0 1 1 1 1 1 1 1 0 0
                            0 0 0 1 1 1 1 1 0 0 0
                            0 0 0 0 1 1 1 0 0 0 0
                            0 0 0 0 0 1 0 0 0 0 0];
                case 2
                    icon = [1 1 1 1 1 1 1 1 1
                            0 1 1 1 1 1 1 1 0
                            0 0 1 1 1 1 1 0 0
                            0 0 0 1 1 1 0 0 0
                            0 0 0 0 1 0 0 0 0];
                case 3
                    icon = [1 1 0 0 0 0 0 0 1 1
                            0 1 1 0 0 0 0 1 1 0
                            0 0 1 1 0 0 1 1 0 0
                            0 0 0 1 1 1 1 0 0 0
                            0 0 0 0 1 1 0 0 0 0
                            0 0 0 1 1 1 1 0 0 0
                            0 0 1 1 0 0 1 1 0 0
                            0 1 1 0 0 0 0 1 1 0
                            1 1 0 0 0 0 0 0 1 1];
                case 4
                    icon = [1 0 0 0 0 0
                            1 1 0 0 0 0
                            1 1 1 0 0 0
                            1 1 1 1 0 0
                            1 1 1 1 1 0
                            1 1 1 1 1 1
                            1 1 1 1 1 0
                            1 1 1 1 0 0
                            1 1 1 0 0 0
                            1 1 0 0 0 0
                            1 0 0 0 0 0];
            end
        end
    end
end