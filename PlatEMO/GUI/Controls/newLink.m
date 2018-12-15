classdef newLink < newGUI
%newLabelLink - Hyperlink label.
%
%   h = newLink(Parent,Pos,Str,...) creates a hyperlink label which has a
%   parent of Parent, a position of Pos and a text of Str. Str stands for
%   the URL which will be opened when click on the label.
%
%   h.callback can be modified by the user, then other function can be
%   executed when clicking on the label, instead of visiting the URL.
%
%   This control is a javax.swing based object, so the properties of it
%   cannot be set as well as uicontrol. Use h.jLabel to get the
%   java.swing.JLabel object.
%
%   Example:
%       newLink(f,[10 10 100 100,1 0 1 0],'www.mathworks.com')
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
        jLabel;     % The text object (java)
        hLabel;     % The text object (matlab)
    end
    methods
        %% Constructor
        function obj = newLink(parent,pos,str,varargin)
            handle = uipanel('Title','','BorderType','none','BackgroundColor',parent.color,...
                             'Parent',parent.handle,'Units','pixels','Position',pos(1:4));
            obj@newGUI(3,parent,handle,pos(5:8),'callback',@(~,~)web(['http://',str],'-browser'),varargin{:});
            % Create the text
            obj.jLabel = javaObjectEDT('javax.swing.JLabel',['<html><u><FONT color="blue">',str,'</FONT></html>']);
            [obj.jLabel,obj.hLabel] = javacomponent(obj.jLabel,[1 1 pos(3:4)],obj.handle);
            color = parent.color;
            obj.jLabel.setBackground(javaObjectEDT('java.awt.Color',color(1),color(2),color(3)));
            obj.jLabel.setFont(javaObjectEDT('java.awt.Font',obj.handle.FontName,0,obj.handle.FontSize));
            % Update the size of obj.hLabel by obj.handle.Position
            obj.newListener(obj.handle,'SizeChanged',@(~,~)set(obj.hLabel,'Position',[1,1,obj.handle.Position(3:4)]));
            % Update obj.jLabel.Enable
            obj.newListener(obj,'state','PostSet',@obj.updateEnable);
        end
    end
    methods(Access = ?newGUI)
        function moveIn(obj)
            obj.figure.handle.Pointer = 'hand';
        end
        function moveOut(obj)
            obj.figure.handle.Pointer = 'arrow';
        end
    end
    methods(Access = protected)
        %% Update obj.jLabel.Enable
        function updateEnable(obj,hObject,eventdata)
            if ~isempty(obj.jLabel)
                obj.jLabel.setEnabled(obj.state)
            end
        end
    end
end