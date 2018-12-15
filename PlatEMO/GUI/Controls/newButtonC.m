classdef newButtonC < newButton
%newButtonC - Combined button.
%
%   h = newButtonC(Parent,Pos,Icon,Str,...) creates a combined button which
%   has a parent of Parent, a position of Pos, an icon of Icon and a text
%   of Str. This type of button has an icon at the top of a text. The text
%   on the button can have multiple lines (set as string cell).
%
%   Example:
%       newButtonC(f,[10 10 100 100,1 0 1 0],cdata,'test')
%
%   See also newGUI, newButton

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetObservable)
        moveoutcallback = @(varargin)[];
    end
    properties(SetAccess = protected)
        icon;           % The icon data
    end
    methods
        %% Constructor
        function obj = newButtonC(parent,pos,icon,str,varargin)
            if ~iscell(str)
                str = {str};
            end
            obj@newButton(parent,pos,['<html><center>',strjoin(str,'<br>'),'</center></html>'],'FontSize',9,varargin{:});
            obj.icon = icon;
            obj.updateCData();
        end
    end
    methods(Access = ?newGUI)
        function moveOut(obj)
            obj.moveoutcallback(obj,[]);
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
                CData(ceil(pos(4)*0.25)+(-sizeI(1)+1:size(obj.icon,1)-sizeI(1)),...
                      ceil(pos(3)*0.50)+(-sizeI(2):size(obj.icon,2)-sizeI(2)-1),:) = obj.icon;
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