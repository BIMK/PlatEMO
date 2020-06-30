classdef newMenu_Item < newButton
%newMenu_Item - The object used in newMenu.
%
%   This is the class of item used in newMenu, which cannot be instantiated
%   independently.
%
%   See also newMenu

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = protected)
        icon;           % The icon data
    end
    methods(Access = ?newMenu)
        %% Constructor
        function obj = newMenu_Item(parent,pos,icon,str,varargin)
            obj@newButton(parent,pos,['<html><tr><td width=9999 align=left><nobr>',repmat('&nbsp',1,5),' ',str],'FontSize',9,varargin{:});
            obj.icon = icon;
            obj.updateCData();
        end
    end
    methods(Access = ?newGUI)
        function mouseUp(obj)
            obj.value = ~obj.value;
            obj.callback(obj,[]);
            obj.parent.visible = false;
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
                      4:3+size(obj.icon,2),:) = obj.icon;
                transparent = repmat(any(isnan(CData),3),1,1,3);
                len         = sum(sum(sum(transparent)))/3;
                CData       = min(max(CData,0),1);
                % The foreground color
                parentColor      = repmat(reshape(obj.parent.color,1,1,3),pos(4),pos(3),1);
                obj.CDataNormal  = CData;
                obj.CDataUnable  = (CData+parentColor)/2;
                obj.CDataMove    = CData;
                obj.CDataPress   = CData;
                obj.CDataUnableP = (CData+parentColor)/2;
                % The background color
                obj.CDataNormal(transparent)  = repmat(obj.parent.color,len,1);
                obj.CDataUnable(transparent)  = repmat(obj.parent.color,len,1);
                obj.CDataMove(transparent)    = repmat((obj.parent.color+[.8 .8 1])/2,len,1);
                obj.CDataPress(transparent)   = repmat((obj.parent.color+[.6 .6 1])/2,len,1);
                obj.CDataUnableP(transparent) = repmat([.9 .9 .9],len,1);
                obj.updateStyle();
            catch
            end
        end
    end
end