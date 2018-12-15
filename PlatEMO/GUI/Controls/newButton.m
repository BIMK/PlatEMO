classdef newButton < newGUI
%newButton - Button.
%
%   h = newButton(Parent,Pos,Str,...) creates a simple button which has a
%   parent of Parent, a position of Pos and a text of Str.
%
%   h.value denotes the times of clicking on the button. If the number of
%   times is even, then h.value = false, otherwise h.value = true. This
%   property is meaningful when h.choosed = true.
%
%   h.choosed = true denotes that it is a toggle button, otherwise it is a
%   general button.
%
%   Example:
%       newButton(f,[10 10 100 100,1 0 1 0],'test')
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

    properties(SetObservable)
        value   = false;   	% Whether the button is being chosen
        choosed = false;    % Whether the button can be chosen
    end
    properties(SetAccess = protected)
        CDataNormal;        % The CData of normal state
        CDataUnable;        % The CData of unable state
        CDataMove;          % The CData of state when mouse is moving in
        CDataPress;         % The CData of state when mouse is pressing
        CDataUnableP;       % The CData of unable state and value = true, choosed = true
    end
    methods
        %% Constructor
        function obj = newButton(parent,pos,str,varargin)
            % Create the button
            handle = uicontrol('ForegroundColor',[.2 .2 .2],'FontName','Microsoft YaHei','FontSize',11,...
                               'Parent',parent.handle,'Style','pushbutton','Enable','inactive','String',str,'Units','pixels','Position',pos(1:4));
            obj@newGUI(3,parent,handle,pos(5:8),varargin{:});
            % Update the button style by obj.moved, obj.pressed, obj.value
            % or obj.state
            obj.newListener(obj,{'moved','pressed','value','state'},'PostSet',@obj.updateStyle);
            % Update the CData by obj.position
            obj.newListener(obj.handle,'SizeChanged',@obj.updateCData);
            obj.updateCData();
        end
    end
    methods(Access = ?newGUI)
        function mouseUp(obj)
            obj.value = ~obj.value;
            obj.callback(obj,[]);
        end
    end
    methods(Access = protected)
        %% Update the CData
        function updateCData(obj,hObject,eventdata)
            try
                pos = floor(obj.handle.Position);
                % The background color
                obj.CDataNormal  = repmat(reshape(obj.parent.color,1,1,3),pos(4),pos(3),1);
                obj.CDataUnable  = obj.CDataNormal;
                obj.CDataMove    = repmat(reshape((obj.parent.color+[.8 .8 1])/2,1,1,3),pos(4),pos(3),1);
                obj.CDataPress   = repmat(reshape((obj.parent.color+[.6 .6 1])/2,1,1,3),pos(4),pos(3),1);
                obj.CDataUnableP = repmat(reshape([.9 .9 .9],1,1,3),pos(4),pos(3),1);
                % The foreground color
                obj.CDataNormal([1,end],:,:) = 0.5;
                obj.CDataNormal(:,[1,end],:) = 0.5;
                obj.CDataUnable([1,end],:,:) = 0.75;
                obj.CDataUnable(:,[1,end],:) = 0.75;
                obj.CDataMove([1,end],:,:)   = 0.5;
                obj.CDataMove(:,[1,end],:)   = 0.5;
                obj.CDataPress([1,end],:,:)  = 0.5;
                obj.CDataPress(:,[1,end],:)  = 0.5;
                obj.updateStyle();
            catch
            end
        end
        %% Update the button style
        function updateStyle(obj,hObject,eventdata)
            if obj.state
                obj.handle.ForegroundColor = [.2 .2 .2];
                if obj.value && obj.choosed || obj.pressed
                    obj.handle.CData = obj.CDataPress;
                elseif obj.moved
                    obj.handle.CData = obj.CDataMove;
                else
                	obj.handle.CData = obj.CDataNormal;
                end
            else
                obj.handle.ForegroundColor = [.5 .5 .5];
                if obj.value && obj.choosed
                    obj.handle.CData = obj.CDataUnableP;
                else
                    obj.handle.CData = obj.CDataUnable;
                end
            end
        end
    end
end