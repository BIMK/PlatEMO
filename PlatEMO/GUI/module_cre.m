classdef module_cre < handle
%module_cre - Creation module.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        GUI;                % The GUI object
        app = struct();     % All the components
        Blocks;             % All blocks
        Graph;              % Adjacent matrix
        Boxes;              % Boxes for all blocks
        Lines;              % Lines for all edges
    end
    methods(Access = ?GUI)
        %% Constructor
        function obj = module_cre(GUI)
            % The main grid
            obj.GUI = GUI;
            obj.app.maingrid = GUI.APP(3,1,uigridlayout(obj.GUI.app.maingrid,'RowHeight',{20,30,'1x'},'ColumnWidth',{'3x',5,1,'1x',1,'1.5x'},'Padding',[5 5 5 5],'RowSpacing',5,'ColumnSpacing',0,'BackgroundColor','w'));
            obj.app.label(1) = GUI.APP(1,1,uilabel(obj.app.maingrid,'Text','Algorithm creation','HorizontalAlignment','center','FontSize',11,'FontColor',[.5 .5 .5]));
            obj.app.label(2) = GUI.APP(1,4,uilabel(obj.app.maingrid,'Text','Problem selection','HorizontalAlignment','center','FontSize',11,'FontColor',[.5 .5 .5]));
            obj.app.label(3) = GUI.APP(1,6,uilabel(obj.app.maingrid,'Text','Training and validation','HorizontalAlignment','center','FontSize',11,'FontColor',[.5 .5 .5]));
            GUI.APP([1 3],3,uipanel(obj.app.maingrid,'BackgroundColor',[.8 .8 .8]));
            GUI.APP([1 3],5,uipanel(obj.app.maingrid,'BackgroundColor',[.8 .8 .8]));

            % The first panel
            obj.app.grid(1)    = GUI.APP(2,1,uigridlayout(obj.app.maingrid,'RowHeight',{1,'1x',1},'ColumnWidth',{18,18,70},'Padding',[5 5 5 5],'RowSpacing',0,'ColumnSpacing',7,'BackgroundColor',[.95 .95 1]));
            tempPanel          = GUI.APP(2,1,uipanel(obj.app.grid(1),'BorderType','none','BackgroundColor',[.95 .95 1]));
            obj.app.buttonA(1) = uibutton(tempPanel,'Position',[-2.5 -2.5 24 24],'Text','','Icon',obj.GUI.icon.loadtable,'BackgroundColor',[.95 .95 1],'Tooltip','Load an algorithm','Interruptible','off','BusyAction','cancel');
            tempPanel          = GUI.APP(2,2,uipanel(obj.app.grid(1),'BorderType','none','BackgroundColor',[.95 .95 1]));
            obj.app.buttonA(2) = uibutton(tempPanel,'Position',[-2.5 -2.5 24 24],'Text','','Icon',obj.GUI.icon.savetable,'BackgroundColor',[.95 .95 1],'Tooltip','Save the problem','Interruptible','off','BusyAction','cancel');
            obj.app.buttonA(3) = GUI.APP([1 3],3,uibutton(obj.app.grid(1),'Text','Add','BackgroundColor',[.95 .95 1],'Tooltip','Check the validity of the algorithm','Interruptible','off','BusyAction','cancel','ButtonpushedFcn',@obj.cb_addblock));
            
            % The second panel
            obj.app.grid(2) = GUI.APP(3,1,uigridlayout(obj.app.maingrid,'RowHeight',{'1x'},'ColumnWidth',{'1x'},'Padding',[0 0 15 10],'BackgroundColor',[.95 .95 1],'RowSpacing',0,'ColumnSpacing',0)); 
            obj.app.canvas  = GUI.APP(1,1,uiaxes(obj.app.grid(2),'NextPlot','add','ClippingStyle','rectangle','XColor','none','YColor','none','view',[0 90],'XLim',[0 100],'YLim',[0 100]));
            
            % The third panel
            
            
            % Initialization
            obj.app.canvas.UserData = struct('type',0,'object',[],'lastpos',[]);
            set(obj.GUI.app.figure,'WindowButtonDownFcn',@obj.cb_figdown,'WindowButtonMotionFcn',@obj.cb_figmove,'WindowButtonUpFcn',@obj.cb_figup);
        end
    end
    methods(Access = private)
        %% Add a new block
        function cb_addblock(obj,~,~)
            obj.Boxes = [obj.Boxes,hggroup(obj.app.canvas,'ButtonDownFcn',@obj.cb_block,'UserData',struct('origin',[],'termi',[]))];
            rectangle(obj.Boxes(end),'Position',[0 90 10 10],'FaceColor',[.95 .95 1],'EdgeColor','k','Curvature',[0.1 0.1],'HitTest',false);
            text(obj.Boxes(end),5,95,['Block ',num2str(length(obj.Boxes))],'HorizontalAlignment','center','Clipping',true,'HitTest',false);
            obj.Graph = [obj.Graph,zeros(size(obj.Graph,1),1);zeros(1,size(obj.Graph,2)+1)];
        end
        %% Press button on the figure
        function cb_figdown(obj,~,~)
            obj.app.canvas.UserData.type = 0;
        end
        %% Move on the figure
        function cb_figmove(obj,~,~)
            if obj.app.canvas.UserData.type == 1
                % Move a block
                h = obj.app.canvas.UserData.object;
                h.Children(2).Position(1:2) = h.Children(2).Position(1:2) + obj.app.canvas.CurrentPoint(1,1:2) - obj.app.canvas.UserData.lastpos;
                h.Children(1).Position(1:2) = h.Children(2).Position(1:2) + 5;
                % Move the curve starting from the block
                for k = h.UserData.origin
                    k.UserData = k.UserData + 0.5*(obj.app.canvas.CurrentPoint(1,1:2)-obj.app.canvas.UserData.lastpos);
                    P = GetCurve(h.Children(1).Position(1:2),[k.Children(2).XData(end),k.Children(2).YData(end)],k.UserData);
                    if P(1,1) < P(end,1) && k.Children(1).String(1) > 999
                        k.Children(1).String = '50%(30) ▶';
                    elseif P(1,1) > P(end,1) && k.Children(1).String(1) <= 999
                        k.Children(1).String = '◀ 50%(30)';
                    end       
                    k.Children(2).XData = P(:,1);
                    k.Children(2).YData = P(:,2);
                    k.Children(1).Rotation      = atan((P(ceil(end/2)+1,2)-P(ceil(end/2),2))./(P(ceil(end/2)+1,1)-P(ceil(end/2),1)))/pi*180;
                    k.Children(1).Position(1:2) = P(ceil(end/2),:);
                end
                % Move the curve ending at the block
                for k = h.UserData.termi
                    k.UserData = k.UserData + 0.5*(obj.app.canvas.CurrentPoint(1,1:2)-obj.app.canvas.UserData.lastpos);
                    P = GetCurve([k.Children(2).XData(1),k.Children(2).YData(1)],h.Children(1).Position(1:2),k.UserData);
                    if P(1,1) < P(end,1) && k.Children(1).String(1) > 999
                        k.Children(1).String = '50%(30) ▶';
                    elseif P(1,1) > P(end,1) && k.Children(1).String(1) <= 999
                        k.Children(1).String = '◀ 50%(30)';
                    end 
                    k.Children(2).XData = P(:,1);
                    k.Children(2).YData = P(:,2);
                    k.Children(1).Rotation      = atan((P(ceil(end/2)+1,2)-P(ceil(end/2),2))./(P(ceil(end/2)+1,1)-P(ceil(end/2),1)))/pi*180;
                    k.Children(1).Position(1:2) = P(ceil(end/2),:);
                end
            elseif obj.app.canvas.UserData.type == 2
                % Move a curve
                h = obj.app.canvas.UserData.object;
                h.UserData = h.UserData + 2*(obj.app.canvas.CurrentPoint(1,1:2)-obj.app.canvas.UserData.lastpos);
                P = GetCurve([h.Children(2).XData(1),h.Children(2).YData(1)],[h.Children(2).XData(end),h.Children(2).YData(end)],h.UserData);
                h.Children(2).XData = P(:,1);
                h.Children(2).YData = P(:,2);
                h.Children(1).Rotation      = atan((P(ceil(end/2)+1,2)-P(ceil(end/2),2))./(P(ceil(end/2)+1,1)-P(ceil(end/2),1)))/pi*180;
                h.Children(1).Position(1:2) = P(ceil(end/2),:);
            end
            obj.app.canvas.UserData.lastpos = obj.app.canvas.CurrentPoint(1,1:2);
        end
        %% Release button on the figure
        function cb_figup(obj,~,~)
            if obj.app.canvas.UserData.type == 0
                ClearHighlight(obj.app.canvas);
                obj.app.canvas.UserData.object = [];
            else
                obj.app.canvas.UserData.type = 0;
            end
        end
        %% Click on a block
        function cb_block(obj,h,event)
            ClearHighlight(obj.app.canvas);
            if ~isempty(obj.app.canvas.UserData.object) && ~isnumeric(obj.app.canvas.UserData.object.UserData) && obj.app.canvas.UserData.object ~= h
                originNo = find(obj.Boxes==obj.app.canvas.UserData.object);
                termiNo  = find(obj.Boxes==h);
                if obj.Graph(originNo,termiNo) == 0
                    % Add a new curve
                    origin    = obj.app.canvas.UserData.object.Children(1).Position;
                    termi     = h.Children(1).Position;
                    obj.Lines = [obj.Lines,hggroup(obj.app.canvas,'ButtonDownFcn',@obj.cb_curve,'UserData',(origin(1:2)+termi(1:2))/2)];
                    P = GetCurve(origin(1:2),termi(1:2),obj.Lines(end).UserData);
                    line(obj.Lines(end),P(:,1),P(:,2),'Color','k');
                    if origin(1) <= termi(1)
                        str = '50%(30) ▶';
                    else
                        str = '◀ 50%(30)';
                    end
                    text(obj.Lines(end),P(ceil(end/2),1),P(ceil(end/2),2),str,'Rotation',atan((P(ceil(end/2)+1,2)-P(ceil(end/2),2))./(P(ceil(end/2)+1,1)-P(ceil(end/2),1)))/pi*180,'Margin',0.1,'BackgroundColor','w','FontSize',10,'HorizontalAlignment','center','Clipping',true,'HitTest',false);
                    obj.app.canvas.Children = obj.app.canvas.Children([2:end,1]);
                    obj.Graph(originNo,termiNo) = 0.5;
                    obj.app.canvas.UserData.object.UserData.origin = [obj.app.canvas.UserData.object.UserData.origin,obj.Lines(end)];
                    h.UserData.termi = [h.UserData.termi,obj.Lines(end)];
                end
            else
                h.Children(2).EdgeColor         = 'b';
                obj.app.canvas.UserData.type    = 1;
                obj.app.canvas.UserData.object  = h;
                obj.app.canvas.UserData.lastpos = obj.app.canvas.CurrentPoint(1,1:2);
            end
        end
        %% Click on a curve
        function cb_curve(obj,h,event)
            ClearHighlight(obj.app.canvas);
            h.Children(2).Color             = 'b';
            obj.app.canvas.UserData.type    = 2;
            obj.app.canvas.UserData.object  = h;
            obj.app.canvas.UserData.lastpos = obj.app.canvas.CurrentPoint(1,1:2);
        end
    end
end

%% Clear the highlighted block or curve
function ClearHighlight(canvas)
    if ~isempty(canvas.UserData.object)
        if isnumeric(canvas.UserData.object.UserData)
            % The current hggroup is a curve
            canvas.UserData.object.Children(2).Color = 'k';
        else
            % The current hggroup is a block
            canvas.UserData.object.Children(2).EdgeColor = 'k';
        end
    end
end

%% Get the points of a Bezier curve
function P = GetCurve(origin,termi,mid)
    t = linspace(0,1,norm(termi-origin)/2)';
    P = (1-t).^2.*origin + 2*t.*(1-t).*mid + t.^2.*termi;
end