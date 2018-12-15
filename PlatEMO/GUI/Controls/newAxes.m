classdef newAxes < newGUI
%newAxes - Axes.
%
%   h = newAxes(Parent,Pos) creates an axes which has a parent of Parent
%   and a position of Pos.
%
%   Example:
%       newAxes(f,[10 10 100 100,1 0 1 0])
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

    properties
        mouseKind = 0;      % The kind of the mouse action
        axisFix   = false;  % Whether the axis range is fixed
        viewFix   = true;   % Whether the view is fixed
    end
    methods
        %% Constructor
        function obj = newAxes(parent,pos)
            handle = axes('NextPlot','add','Box','on','FontName','Times New Roman','FontSize',13,'Parent',parent.handle,'Units','pixels','Position',pos(1:4));
            obj@newGUI(3,parent,handle,pos(5:8));
        end
        %% Plot data (replace)
        function draw(obj,Data,varargin)
            [N,M] = size(Data);
            % The size of the figure
            if obj.handle.Position <= [inf inf 400 300]
                Size = [3 5 .8 8];
            else
                Size = [6 8 2 13];
            end
            % The style of the figure
            if nargin < 3
                if M == 2
                    varargin = {'ok','MarkerSize',Size(1),'Marker','o','Markerfacecolor',[.7 .7 .7],'Markeredgecolor',[.4 .4 .4]};
                elseif M == 3
                    varargin = {'ok','MarkerSize',Size(2),'Marker','o','Markerfacecolor',[.7 .7 .7],'Markeredgecolor',[.4 .4 .4]};
                elseif M > 3
                    varargin = {'Color',[.5 .5 .5],'LineWidth',Size(3)};
                end
            end
            % Draw the figure
            subplot(obj.handle); cla;
            set(obj.handle,'FontSize',Size(4));
            if M == 2
                plot(Data(:,1),Data(:,2),varargin{:});
            elseif M == 3
                plot3(Data(:,1),Data(:,2),Data(:,3),varargin{:});
            elseif M > 3
                Label = repmat([0.99,2:M-1,M+0.01],N,1);
                Data(2:2:end,:)  = fliplr(Data(2:2:end,:));
                Label(2:2:end,:) = fliplr(Label(2:2:end,:));
                Data  = Data';
                Label = Label';
                plot(Label(:),Data(:),varargin{:});
            end
            if ~obj.axisFix
                axis tight;
                if M > 3
                    set(obj.handle,'XLim',[1,M]);
                end
            end
            if ~obj.viewFix
                if M == 2
                    view(0,90);
                    xlabel('\itf\rm_1'); ylabel('\itf\rm_2');
                elseif M == 3
                    view(135,30);
                    xlabel('\itf\rm_1'); ylabel('\itf\rm_2'); zlabel('\itf\rm_3');
                elseif M > 3
                    view(0,90);
                    xlabel('Dimension No.'); ylabel('Value');
                end
            end
        end
    end
    methods(Access = ?newGUI)
        function moveIn(obj)
            switch obj.mouseKind
                case 1
                    zoom(obj.figure.handle,'inmode');
                case 2
                    zoom(obj.figure.handle,'outmode');
                case 3
                    pan(obj.figure.handle,'onkeepstyle');
                case 4
                    rotate3d(obj.figure.handle,'on');
            end
        end
        function moveOut(obj)
            if obj.mouseKind
                % Call the WindowButtonUpFcn of gcf forcibly
                try
                    buttonUpFcn = obj.figure.handle.WindowButtonUpFcn;
                    buttonUpFcn{1}(obj.figure.handle,[],buttonUpFcn{2:end});
                catch
                end
                % Reset the mouse action
                zoom(obj.figure.handle,'off');
                pan(obj.figure.handle,'off');
                rotate3d(obj.figure.handle,'off');
            end
        end
    end
end