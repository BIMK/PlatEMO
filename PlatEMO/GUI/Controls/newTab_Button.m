classdef newTab_Button < newButton
%newTab_Button - The object used in newTab.
%
%   This is the class of toggle button used in newTab, which cannot be
%   instantiated independently.
%
%   See also newTab

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods(Access = ?newTab)
        %% Constructor
        function obj = newTab_Button(parent,pos,str,varargin)
            obj@newButton(parent,pos,str,varargin{:});
            obj.updateCData();
        end
    end
    methods(Access = protected)
        %% Update the CData
        function updateCData(obj,hObject,eventdata)
            try
                pos = floor(obj.handle.Position);
                obj.CDataNormal = repmat(reshape(obj.parent.color,1,1,3),pos(4),pos(3),1);
                obj.CDataPress  = obj.CDataNormal;
                obj.CDataNormal(end,:,:)    = 0.5;
                obj.CDataPress(1,:,:)       = 0.5;
                obj.CDataPress(:,[1,end],:) = 0.5;
                obj.updateStyle();
            catch
            end
        end
        %% Update the button style
        function updateStyle(obj,hObject,eventdata)
            if obj.state
                if obj.value
                    obj.handle.CData = obj.CDataPress;
                    obj.handle.ForegroundColor = [.5 .5 1];
                elseif obj.moved || obj.pressed
                    obj.handle.CData = obj.CDataNormal;
                    obj.handle.ForegroundColor = [.5 .5 1];
                else
                    obj.handle.CData = obj.CDataNormal;
                    obj.handle.ForegroundColor = [.2 .2 .2];
                end
            else
                if obj.value
                    obj.handle.CData = obj.CDataPress;
                    obj.handle.ForegroundColor = [.8 .8 .8];
                else
                    obj.handle.CData = obj.CDataNormal;
                    obj.handle.ForegroundColor = [.8 .8 .8];
                end
            end
        end
    end
end