classdef module < handle
%module - The superclass of all modules.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = protected, SetObservable)
        GUI;                    % The main GUI object
        panel;                  % The main panel of this module
        control = struct();     % All the controls
        data    = [];           % For storing the data
    end
    methods(Access = protected)
        %% Create the module object
        function obj = module(GUI,panel)
            obj.GUI   = GUI;
            obj.panel = panel;
        end
    end
    methods(Static)
        %% Obtain the summary
        function words = summary()
            words = [];
        end
        %% Callback of the figure in guidance mode
        function guidance(obj)
            obj.GUI.figure.busy   = false;
            obj.GUI.figure.enable = true;
            delete(obj.control.guideLabel);
        end
    end
    methods
        %% Show the guidance
        function showGuidance(obj)
            obj.control.guideLabel = newLabel(obj.GUI.figure,[1 1 1 20,1 0 1 0],'','FontSize',12,'ForegroundColor',[1 1 1],'BackgroundColor',[1 .1 .1],'HorizontalAlignment','left');
            obj.control.guideIndex = 0;
            obj.GUI.figure.handle.WindowButtonUpFcn = @(~,~)obj.guidance(obj);
            obj.guidance(obj);
        end
        function updateGuideLabel(obj,str,target)
            target.state = true;
            obj.control.guideLabel.handle.String = str;
            if target.absPos(1)+obj.control.guideLabel.handle.Extent(3) <= obj.GUI.figure.position(3)
                obj.control.guideLabel.position(1:3) = [target.absPos(1),target.absPos(2)+target.absPos(4)+15,obj.control.guideLabel.handle.Extent(3)];
            else
                obj.control.guideLabel.position(1:3) = [target.absPos(1)+target.absPos(3)-obj.control.guideLabel.handle.Extent(3),target.absPos(2)+target.absPos(4)+15,obj.control.guideLabel.handle.Extent(3)];
            end
        end
        %% Calculate the specified metric value of a population
        function value = Metric(obj,metric,PopObj,PF)
            if isa(PopObj,'INDIVIDUAL')
                PopObj = PopObj(all(PopObj.cons<=0,2)).objs;
            end
            NonDominated = NDSort(PopObj,1) == 1;
            try
                value = metric(PopObj(NonDominated,:),PF);
            catch
                value = NaN;
            end
        end
    end
end