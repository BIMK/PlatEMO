classdef newMpanel < newPanel
%newMpanel - Set of multiple panels.
%
%   h = newMpanel(Parent,Pos,N) creates a panel which has a parent of
%   Parent and a position of Pos. The panel contains N subpanels.
%
%   Use h.panels(i) to obtain the i-th subpanel in h.
%
%   Example:
%       newMpanel(f,[10 10 100 100,1 0 1 0],3)
%
%   See also newGUI, newPanel

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = protected)
        panels;     % All the sub panels
    end
    methods
        %% Constructor
        function obj = newMpanel(parent,pos,nPanel)
            obj@newPanel(parent,pos,parent.color);
            % Create the sub panels
            for i = 1 : nPanel
                obj.panels = [obj.panels,newPanel(obj,[(i-1)*pos(3)/nPanel,0,pos(3)/nPanel,pos(4),0 0 1 1],parent.color)];
                if i > 1
                    newLine(obj.panels(end),[1,6,1,pos(4)-8,1 0 1 1]);
                end
            end
        end
        %% Move the gap
        function move(obj,index,offset)
            if obj.panels(index).position(3) + offset(1) > 50 && obj.panels(index+1).position(3) - offset(1) > 50
                obj.panels(index).position(3)       = obj.panels(index).position(3) + offset(1);
                obj.panels(index+1).position([1,3]) = [obj.panels(index+1).position(1)+offset(1),obj.panels(index+1).position(3)-offset(1)];
            end
        end
    end
end