classdef newEditPass < newGUI
%newEditPass - Password edit box.
%
%	h = newEditPass(Parent,Pos,...) creates a password edit box which
%   has a parent of Parent and a position of Pos.
%
%   Use h.string to get the text of the edit box.
%
%   This control is a javax.swing based object, so the properties of it
%   cannot be set as well as uicontrol. Use h.jEdit to get the
%   corresponding javax.swing.JPasswordField object.
%
%	Example:
%       newEditPass(f,[10 10 100 100,1 0 1 0])
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

    properties(SetAccess = protected)
        jEdit;      % The edit box object (java)
        hEdit;      % The edit box object (matlab)
    end
    methods
        %% Constructor
        function obj = newEditPass(parent,pos,varargin)
            handle = uipanel('Title','','BorderType','none','BackgroundColor',parent.color,...
                             'Parent',parent.handle,'Units','pixels','Position',pos(1:4));
            obj@newGUI(3,parent,handle,pos(5:8),varargin{:});
            % Create the edit box
            obj.jEdit = javaObjectEDT('javax.swing.JPasswordField');
            [obj.jEdit,obj.hEdit] = javacomponent(obj.jEdit,[1 1 pos(3:4)],obj.handle);
            obj.jEdit.setFont(javaObjectEDT('java.awt.Font',obj.handle.FontName,0,obj.handle.FontSize));
            % Update the size of obj.hEdit by obj.handle.Position
            obj.newListener(obj.handle,'SizeChanged',@(~,~)set(obj.hEdit,'Position',[1,1,obj.handle.Position(3:4)]));
            % Update obj.jEdit.Enable
            obj.newListener(obj,'state','PostSet',@obj.updateEnable);
        end
        %% Get the real string of the object
        function str = string(obj)
            str = obj.jEdit.getPassword()';
        end
    end
    methods(Access = protected)
        %% Update obj.jEdit.Enable
        function updateEnable(obj,hObject,eventdata)
            if ~isempty(obj.jEdit)
                obj.jEdit.setEnabled(obj.state)
            end
        end
    end
end