classdef newBar < newPanel
%newBar - Progress bar.
%
%   h = newBar(Parent,Pos) creates a horizontal progress bar which has a
%   parent of Parent and a position of Pos.
%
%   h.value denotes the current value of the bar, which has a minimum value
%   of 0 and a maximum value of 1.
%
%   Example:
%       newBar(f,[10 10 100 100,1 0 1 0])
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

    properties(SetObservable)
        value = 0;      % The current valu of the Pbar
    end
    properties(SetAccess = protected)
        blocks;         % The blocks in the Pbar
    end
    methods
        %% Constructor
        function obj = newBar(parent,pos)
            pos(3) = round((pos(3)-6)/17)*17+6;
            obj@newPanel(parent,pos,[]);
            for i = 1 : (pos(3)-6)/17
                obj.blocks = [obj.blocks,uicontrol('Parent',obj.handle,'Style','text','Units','pixels','Position',[5+(i-1)*17,4,13,pos(4)-8],'Visible','off','BackgroundColor',[0.078 0.169 0.549])];
            end
            % Update the position of the blocks
            obj.newListener(obj.handle,'SizeChanged',@obj.updateBlock);
        end
        %% The set method of obj.value
        function set.value(obj,value)
            obj.value = min(max(value,0),1);
            % Update the visibilities of the labels by obj.value
            nShown = round(obj.value*length(obj.blocks));
            if nShown > 0 
                [obj.blocks(1:nShown).Visible] = deal('on');
            end
            if nShown < length(obj.blocks)
                [obj.blocks(nShown+1:end).Visible] = deal('off');
            end
        end
    end
    methods(Access = protected)
        function updateBlock(obj,hObject,eventdata)
            if length(obj.blocks) > 1
                Position = cell2mat(get(obj.blocks,'Position'));
            else
                Position = obj.blocks.Position;
            end
            Position(:,[1,3]) = Position(:,[1,3]).*obj.position(3)./obj.oldAbsPos(3);
            Position(:,4)     = max(1,Position(:,4)+obj.position(4)-obj.oldAbsPos(4));
            set(obj.blocks,{'Position'},num2cell(Position,2));
        end
    end
end