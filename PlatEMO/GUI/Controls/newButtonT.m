classdef newButtonT < newButton
%newButton - Tool button.
%
%   h = newButtonT(Parent,Pos,Input,Tip...) creates a button which has a
%   parent of Parent, a position of Pos and a tooltip text of Tip. If Input
%   is a string, then the button has a text; if Input is an image, then the
%   button has an icon.
%
%   Example:
%       newButtonT(f,[10 10 100 100,1 0 1 0],'test')
%       newButtonT(f,[10 10 100 100,1 0 1 0],cdata)
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

    properties(SetAccess = protected)
        icon;           % The icon data
        tip;            % The tip label
    end
    methods
        %% Constructor
        function obj = newButtonT(parent,pos,input,tip,varargin)
            if ischar(input)
                str  = input;
                icon = [];
            else
                str  = '';
                icon = input;
            end
            obj@newButton(parent,pos,str,varargin{:});
            obj.icon = icon;
            tiplabel = uicontrol('Parent',obj.figure.handle,'Style','text','Units','pixels','Position',[0 0 1 1],'FontSize',10,'String',tip,'BackgroundColor',[1 1 .9]);
            obj.tip  = uipanel('Parent',obj.figure.handle,'Units','pixels','Position',[tiplabel.Extent(1:3),16],'HighlightColor',[.5 .5 .5],'BorderType','line','Visible','off');
            set(tiplabel,'Parent',obj.tip,'Position',[0 3 tiplabel.Extent(3) 13]);
            obj.updateCData();
        end
    end
    methods(Access = ?newGUI)
        function moveIn(obj)
            obj.tip.Position(1:2) = [obj.absPos(1)+15,obj.absPos(2)-20];
            obj.tip.Visible = 'on';
        end
        function moveOut(obj)
            obj.tip.Visible = 'off';
        end
    end
    methods(Access = protected)
        %% Update the CData
        function updateCData(obj,hObject,eventdata)
            try
                if isempty(obj.icon)
                    pos = floor(obj.handle.Position)+2;
                    % The background color
                    obj.CDataNormal  = repmat(reshape(obj.parent.color,1,1,3),pos(4),pos(3),1);
                    obj.CDataUnable  = obj.CDataNormal;
                    obj.CDataMove    = repmat(reshape((obj.parent.color+[.8 .8 1])/2,1,1,3),pos(4),pos(3),1);
                    obj.CDataPress   = repmat(reshape((obj.parent.color+[.6 .6 1])/2,1,1,3),pos(4),pos(3),1);
                    obj.CDataUnableP = repmat(reshape([.9 .9 .9],1,1,3),pos(4),pos(3),1);
                    obj.updateStyle();
                else
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
                end
            catch
            end
        end
    end
end