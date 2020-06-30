classdef newTable < newGUI
%newTable - Table.
%
%   h = newTable(Parent,Pos,...) creates a table which has a parent of
%   Parent and a position of Pos.
%
%   Example:
%       newTable(f,[10 10 100 100,1 0 1 0])
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
        jTable;      % The table object (java)
    end
    methods
        %% Constructor
        function obj = newTable(parent,pos,varargin)
            % Create the table
            handle = uitable('Parent',parent.handle,'Units','pixels','Position',pos(1:4));
            obj@newGUI(3,parent,handle,pos(5:8),varargin{:});
        end
        %% Set or get the selected cell
        function [row,column] = selected(obj,row,column)
            if isempty(obj.jTable)
                obj.jTable = obj.findJTable();
            end
            if nargin > 2
                drawnow();
                try
                    obj.jTable.changeSelection(row-1,column-1,false,false);
                catch
                end
            else
                row    = obj.jTable.getSelectedRow + 1;
                column = obj.jTable.getSelectedColumn + 1;
            end
        end
    end
    methods(Access = protected)
        %% Get the corresponding java object
        function jTable = findJTable(obj)
            warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
            obj.handle.TooltipString = 'newTableIdentifier';
            drawnow();
            jTable = obj.findJTable2(obj.figure.handle.JavaFrame.getFigurePanelContainer.getComponent(0));
            obj.handle.TooltipString = '';
        end
        function jTable = findJTable2(obj,jParent)
            try
                jTable = [];
                if strcmp(char(jParent.getToolTipText),'newTableIdentifier')
                    jTable = jParent.getParent.getView.getParent.getParent.getViewport.getView;
                else
                    for i = 0 : jParent.getComponentCount-1
                        jTable = obj.findJTable2(jParent.getComponent(i));
                        if ~isempty(jTable)
                            break;
                        end
                    end
                end
            catch
            end
        end
    end
end