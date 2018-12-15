classdef newLine < newGUI
%newLine - Line.
%
%   h = newLine(Parent,Pos) creates a line which has a parent of Parent and
%   a position of Pos. The height or width should be 1.
%
%   Example:
%       newLine(f,[10 10 100 1,1 0 1 0])
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
        function obj = newLine(parent,pos)
            handle = uipanel(parent.handle,'Title','','HighlightColor',[.5 .5 .5],'BorderType','line','Units','pixels','Position',pos(1:4));
            obj@newGUI(3,parent,handle,pos(5:8));
        end
    end
end