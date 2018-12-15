classdef newPanel < newGUI
%newPanel - Panel.
%
%   h = newPanel(Parent,Pos,Color,...) creates a panel which has a parent
%   of Parent, a position of Pos and a background color of Color. If Color
%   is empty, the panel will have a frame and be transparent.
%
%   Example:
%       newPanel(f,[10 10 100 100,1 0 1 0],[])
%       newPanel(f,[10 10 100 100,1 0 1 0],rand(1,3))
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
        function obj = newPanel(parent,pos,color,varargin)
            if isempty(color)
                handle = uipanel('Title','','HighlightColor',[.5 .5 .5],'BorderType','line','BackgroundColor',parent.color,...
                                 'Parent',parent.handle,'Units','pixels','Position',pos(1:4));
                absOffset = 0;
            else
                handle = uipanel('Title','','BorderType','none','BackgroundColor',color,...
                             'Parent',parent.handle,'Units','pixels','Position',pos(1:4));
                absOffset = -1;
            end
            obj@newGUI(2,parent,handle,pos(5:8),varargin{:},'absOffset',absOffset);
        end
    end
end