classdef module_cre < handle
%module_cre - Creation module.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        GUI;                    % The GUI object
        app       = struct();	% All the components
        blockList;              % List of candidate blocks
        Graph     = digraph();	% All blocks and lines
        data      = {};     	% All the results
        dataTrain = {};         % All the results of training
    end
    methods(Access = ?GUI)
        %% Constructor
        function obj = module_cre(GUI)
            % The main grid
            obj.GUI = GUI;
            obj.app.maingrid = GUI.APP(3,1,uigridlayout(obj.GUI.app.maingrid,'RowHeight',{20,30,'1.2x',1,15,'1x'},'ColumnWidth',{'3x',5,1,'1.1x',1,5,'1.5x'},'Padding',[5 5 5 5],'RowSpacing',5,'ColumnSpacing',0,'BackgroundColor','w'));
            obj.app.label(1) = GUI.APP(1,1,uilabel(obj.app.maingrid,'Text','Algorithm creation','HorizontalAlignment','center','FontSize',11,'FontColor',[.5 .5 .5]));
            obj.app.label(2) = GUI.APP(1,4,uilabel(obj.app.maingrid,'Text','Problem selection','HorizontalAlignment','center','FontSize',11,'FontColor',[.5 .5 .5]));
            obj.app.label(3) = GUI.APP(1,7,uilabel(obj.app.maingrid,'Text','Training','HorizontalAlignment','center','FontSize',11,'FontColor',[.5 .5 .5]));
            obj.app.label(4) = GUI.APP(5,7,uilabel(obj.app.maingrid,'Text','Testing','HorizontalAlignment','center','FontSize',11,'FontColor',[.5 .5 .5]));
            GUI.APP([1 6],3,uipanel(obj.app.maingrid,'BackgroundColor',[.8 .8 .8]));
            GUI.APP([1 6],5,uipanel(obj.app.maingrid,'BackgroundColor',[.8 .8 .8]));
            GUI.APP(4,7,uipanel(obj.app.maingrid,'BackgroundColor',[.8 .8 .8]));

            % The first panel
            obj.app.loadmenu = uicontext(obj.GUI.app.figure,115);
            obj.app.loadmenu.add('  Load a block','',@obj.cb_loadblock);
            obj.app.loadmenu.add('  Load an algorithm','',@obj.cb_loadalgorithm);
            obj.app.loadmenu.flush();
            obj.app.grid(1)    = GUI.APP(2,1,uigridlayout(obj.app.maingrid,'RowHeight',{1,'1x',1},'ColumnWidth',{18,24,18,18,1,18,18,'1x',70,'0.9x'},'Padding',[10 5 10 5],'RowSpacing',0,'ColumnSpacing',7,'BackgroundColor',[.95 .95 1]));
            tempPanel          = GUI.APP(2,1,uipanel(obj.app.grid(1),'BorderType','none','BackgroundColor',[.95 .95 1]));
            obj.app.buttonA(1) = uibutton(tempPanel,'Position',[-2.5 -2.5 24 24],'Text','','Icon',obj.GUI.icon.block,'BackgroundColor',[.95 .95 1],'Tooltip','Add a block','Interruptible','off','BusyAction','cancel','ButtonpushedFcn',@obj.cb_addblock);
            tempPanel          = GUI.APP(2,2,uipanel(obj.app.grid(1),'BorderType','none','BackgroundColor',[.95 .95 1]));
            obj.app.buttonA(2) = uibutton(tempPanel,'Position',[-2.5 -2.5 31 24],'Text','','Icon',obj.GUI.icon.loadtable2,'BackgroundColor',[.95 .95 1],'Tooltip','Load a block or algorithm','Interruptible','off','BusyAction','cancel','ButtonpushedFcn',@(~,~)obj.app.loadmenu.show());
            tempPanel          = GUI.APP(2,3,uipanel(obj.app.grid(1),'BorderType','none','BackgroundColor',[.95 .95 1]));
            obj.app.buttonA(3) = uibutton(tempPanel,'Position',[-2.5 -2.5 24 24],'Text','','Icon',obj.GUI.icon.savetable,'BackgroundColor',[.95 .95 1],'Tooltip','Save the algorithm','Interruptible','off','BusyAction','cancel','ButtonpushedFcn',@obj.cb_savealgorithm);
            tempPanel          = GUI.APP(2,4,uipanel(obj.app.grid(1),'BorderType','none','BackgroundColor',[.95 .95 1]));
            obj.app.buttonA(4) = uibutton(tempPanel,'Position',[-2.5 -2.5 24 24],'Text','','Icon',obj.GUI.icon.code,'BackgroundColor',[.95 .95 1],'Tooltip','Generate code','Interruptible','off','BusyAction','cancel','ButtonpushedFcn',@obj.cb_generatecode);
            GUI.APP([1 3],5,uipanel(obj.app.grid(1),'BackgroundColor',[.8 .8 .8]));
            tempPanel          = GUI.APP(2,6,uipanel(obj.app.grid(1),'BorderType','none','BackgroundColor',[.95 .95 1]));
            obj.app.buttonA(5) = uibutton(tempPanel,'Position',[-2.5 -2.5 24 24],'Text','','Icon',obj.GUI.icon.feature,'BackgroundColor',[.95 .95 1],'Tooltip','Select features to display','Interruptible','off','BusyAction','cancel','ButtonpushedFcn',@obj.cb_selectfeature);
            tempPanel          = GUI.APP(2,7,uipanel(obj.app.grid(1),'BorderType','none','BackgroundColor',[.95 .95 1]));
            obj.app.buttonA(6) = uibutton(tempPanel,'Position',[-2.5 -2.5 24 24],'Text','','Icon',obj.GUI.icon.arrange,'BackgroundColor',[.95 .95 1],'Tooltip','Arrange the blocks','Interruptible','off','BusyAction','cancel','ButtonpushedFcn',@obj.cb_arrange);
            obj.app.buttonA(7) = GUI.APP([1 3],9,uibutton(obj.app.grid(1),'Text','Validation','BackgroundColor',[.95 .95 1],'Tooltip','Check the validity of the algorithm','Interruptible','off','BusyAction','cancel','ButtonpushedFcn',@obj.cb_validation));
            obj.app.dropA      = GUI.APP([1 3],10,uidropdown(obj.app.grid(1),'BackgroundColor',[.95 .95 1],'Tooltip','Select predefined examples','Items',{'User-defined algorithm','General algorithm','Coevolutionary algorithm','Complex algorithm'},'ItemsData',1:4,'Value',1,'Interruptible','off','BusyAction','cancel','ValueChangedFcn',@obj.cb_selectAlgorithm));
            
            % The second panel
            obj.app.grid(2) = GUI.APP([3 6],1,uigridlayout(obj.app.maingrid,'RowHeight',{'1x'},'ColumnWidth',{'1x'},'Padding',[0 0 15 10],'BackgroundColor',[.95 .95 1],'RowSpacing',0,'ColumnSpacing',0)); 
            obj.app.canvas  = GUI.APP(1,1,uiaxes(obj.app.grid(2),'NextPlot','add','ClippingStyle','rectangle','XColor','none','YColor','none','view',[0 90],'XLim',[0 100],'YLim',[0 100],'ButtonDownFcn',@(h,~)HighlightObject(h)));
            obj.app.menu    = uicontext(obj.GUI.app.figure,85);
            obj.app.menu.add('Details',obj.GUI.icon.feature,@obj.cb_checkdetail);
            obj.app.menu.add('Color',obj.GUI.icon.color,@obj.cb_setcolor);
            obj.app.menu.add('Open file',obj.GUI.icon.file,@obj.cb_openfile);
            obj.app.menu.add('Save block',obj.GUI.icon.savetable,@obj.cb_saveblock);
            obj.app.menu.add('Delete',obj.GUI.icon.delete,@obj.cb_delete);
            set(obj.app.menu.gaps([2,4]),'Visible',false);
            obj.app.menu.flush();
            
            % The third panel
            obj.app.grid(3)   = GUI.APP([2 6],4,uigridlayout(obj.app.maingrid,'RowHeight',{16,21,16,21,21,16,21,21,21,4,18,'1.1x','1x'},'ColumnWidth',{'1x','1.1x','1x'},'Padding',[8 10 8 0],'RowSpacing',3,'ColumnSpacing',5,'BackgroundColor','w'));
            [obj.app.stateC,obj.app.labelC] = GUI.GenerateLabelButton(obj.app.grid(3),[1 0 0 1,zeros(1,13)],@obj.cb_filter);
            obj.app.labelC(4) = GUI.APP(11,[1 2],uilabel(obj.app.grid(3),'Text','Problems','FontSize',13,'FontColor',[.9 .5 .2],'FontWeight','bold'));
            obj.app.labelC(5) = GUI.APP(11,3,uilabel(obj.app.grid(3),'HorizontalAlignment','right','FontSize',10,'FontColor',[.9 .5 .2]));
            obj.app.listC     = GUI.APP(12,[1 3],uilistbox(obj.app.grid(3),'FontColor',[.9 .5 .2]));
            obj.app.dropC     = GUI.APP(11,2,uidropdown(obj.app.grid(3),'BackgroundColor','w','FontColor',[.9 .5 .2],'Items',{'All year'},'ValueChangedFcn',@(h,~)GUI.UpdateAlgProListYear(obj.app.listC,h,obj.app.labelC(5),obj.GUI.proList)));
            obj.app.listD     = uilist(obj.app.grid(3),obj.GUI.app.figure,obj.GUI.icon);
            obj.app.listC.ValueChangedFcn = @(~,~)GUI.UpdateAlgProPara(obj.GUI.app.figure,obj.app.listC,obj.app.listD,'PROBLEM',-1);
            obj.app.listD.grid.Padding    = [0,0,0,0];
            GUI.APP(13,[1 3],obj.app.listD.grid);
            
            % The fourth panel
            obj.app.grid(4)    = GUI.APP([2 3],7,uigridlayout(obj.app.maingrid,'RowHeight',{22,22,22,'1x',27},'ColumnWidth',{'1.3x','0.7x','1.3x',5,'1.3x','0.7x','1.3x'},'Padding',[5 5 5 0],'BackgroundColor','w','RowSpacing',5,'ColumnSpacing',5)); 
            GUI.APP(1,[1 2],uilabel(obj.app.grid(4),'Text','Population size','Tooltip','Population size of the trainer'));
            obj.app.editD(1)   = GUI.APP(1,3,uieditfield(obj.app.grid(4),'numeric','Value',50,'Limits',[10 inf],'RoundFractionalValues',true,'Tooltip','Population size of the trainer'));
            GUI.APP(1,[5 6],uilabel(obj.app.grid(4),'Text','Max evaluations','Tooltip','Maximum number of function evaluations of the trainer'));
            obj.app.editD(2)   = GUI.APP(1,7,uieditfield(obj.app.grid(4),'numeric','Value',10000,'ValueDisplayFormat','%d','Limits',[100 inf],'RoundFractionalValues',true,'Tooltip','Maximum number of function evaluations of the trainer'));
            GUI.APP(2,[1 2],uilabel(obj.app.grid(4),'Text','Execution times','Tooltip','Times of using each algorithm to solve the problem in algorithm performance evaluation'));
            obj.app.editD(3)   = GUI.APP(2,3,uieditfield(obj.app.grid(4),'numeric','Value',3,'Limits',[1 inf],'RoundFractionalValues',true,'Tooltip','Times of using each algorithm to solve the problem in algorithm performance evaluation'));
            GUI.APP(2,[5 6],uilabel(obj.app.grid(4),'Text','Parallelization','Tooltip','Perform the training with multiple CPUs'));
            obj.app.checkD     = GUI.APP(2,7,uicheckbox(obj.app.grid(4),'Text','','Tooltip','Perform the training with multiple CPUs','Enable',~isempty(ver('parallel'))));
            GUI.APP(3,1,uilabel(obj.app.grid(4),'Text','File path','Tooltip','The population will be automatically loaded from file before training and saved to file after each iteration'));
            obj.app.editD(4)   = GUI.APP(3,[3 7],uieditfield(obj.app.grid(4),'Value',fullfile(cd,'Algorithms','Blocks','myAlgorithm.mat'),'Tooltip','The population will be automatically loaded from file before training and saved to file after each iteration'));
            obj.app.buttonD(1) = GUI.APP(3,2,uibutton(obj.app.grid(4),'Text','...','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_filepath,obj.app.editD(4)},'Tooltip','The population will be automatically loaded from file before training and saved to file after each iteration'));
            obj.app.axesD      = GUI.APP(4,[1 7],uiaxes(obj.app.grid(4),'BackgroundColor','w','Box','on','FontName','Times New Roman','FontSize',11));
            axtoolbar(obj.app.axesD);
            obj.app.buttonD(2) = GUI.APP(5,3,uibutton(obj.app.grid(4),'push','Text','Start','ButtonpushedFcn',@obj.cb_train));
            obj.app.buttonD(3) = GUI.APP(5,5,uibutton(obj.app.grid(4),'push','Text','Stop','Enable',false,'ButtonpushedFcn',@(~,~)set(obj.app.buttonD(2:3),{'Enable','Text'},{true,'Start';false,'Stop'})));
            obj.app.labelD     = GUI.APP(5,[1 2],uilabel(obj.app.grid(4),'Text','0.00%','WordWrap',true));
            
            % The fifth panel
            obj.app.grid(5)    = GUI.APP(6,7,uigridlayout(obj.app.maingrid,'RowHeight',{'1x',27},'ColumnWidth',{'1.3x','0.7x','1.3x',5,'1.3x','0.7x','1.3x'},'Padding',[5 5 5 0],'BackgroundColor','w','RowSpacing',5,'ColumnSpacing',5)); 
            obj.app.axesE      = GUI.APP(1,[1 7],uiaxes(obj.app.grid(5),'BackgroundColor','w','Box','on'));
            axtoolbar(obj.app.axesE);
            obj.app.buttonE(1) = GUI.APP(2,3,uibutton(obj.app.grid(5),'push','Text','Start','ButtonpushedFcn',@obj.cb_test));
            obj.app.buttonE(2) = GUI.APP(2,5,uibutton(obj.app.grid(5),'push','Text','Stop','Enable',false,'ButtonpushedFcn',@(~,~)set(obj.app.buttonE(1:2),{'Enable','Text'},{true,'Start';false,'Stop'})));
            obj.app.labelE     = GUI.APP(2,[1 2],uilabel(obj.app.grid(5),'Text','0.00%','WordWrap',true));
            obj.app.menuE      = uicontext(obj.GUI.app.figure,120);
            obj.app.menuE.add('  Save best solutions','',@(~,~)GUI.SavePopulation(obj.GUI.app.figure,obj.data{1}.result{end},1));
            obj.app.menuE.add('  Save all solutions','',@(~,~)GUI.SavePopulation(obj.GUI.app.figure,obj.data{1}.result{end},2));
            obj.app.menuE.flush();
            obj.app.buttonE(3) = GUI.APP(2,7,uibutton(obj.app.grid(5),'push','Text','Save','Enable',false,'ButtonpushedFcn',@(~,~)obj.app.menuE.show()));
            
            % Initialization
            obj.app.canvas.UserData = struct('features',[0 1 0 0 0 0 1 0],'object',[],'lastpos',[],'connectable',false);
            
            obj.cb_filter([],[],0);
            obj.app.listC.Value = 'SOP_F1';
            GUI.UpdateAlgProPara(obj.GUI.app.figure,obj.app.listC,obj.app.listD,'PROBLEM',-1);
            obj.app.listD.items(1).edit(3).Value = '10';
            
            % Read block list
            obj.blockList = {};
            Folders = split(genpath(fullfile(fileparts(mfilename('fullpath')),'..','Algorithms','Blocks')),pathsep);
            for i = 1 : length(Folders) - 1
                Files = what(Folders{i});
                Files = Files.m;
                for j = 1 : length(Files)
                    try
                        f = fopen(Files{j});
                        if contains(fgetl(f),'< BLOCK')
                            Comment   = {};
                            Parameter = {};
                            str       = fgetl(f);
                            while ischar(str) && ~isempty(regexp(str,'^\s*%\s*','once'))
                                str = regexprep(str,'^\s*%\s*','','once');
                                str = regexp(str,'\s*---\s*','split');
                                if length(str) < 2
                                    Comment = [Comment,str];
                                else
                                    if length(str) == 2
                                        str{3} = '';
                                    end
                                    Parameter = [Parameter;str(1:3)];
                                end
                                str = fgetl(f);
                            end
                            obj.blockList = [obj.blockList;Files{j}(1:end-2),strjoin(Comment,' '),{Parameter}];
                        end
                    catch
                    end
                end
            end
        end
    end
    methods(Access = private)
        %% Add a new block
        function cb_addblock(obj,~,~)
            HighlightObject(obj.app.canvas);
            P = obj.GUI.app.figure.Position;
            C = obj.GUI.app.figure.CurrentPoint;
            f = uifigure('Name','Add a new block','Position',[P(1)+C(1) P(2)+C(2)-200 300 180],'icon',obj.GUI.icon.logo1,'WindowStyle','modal','Resize','off');
            temp.tip    = uilabel(f,'Position',[160 140 120 40],'FontSize',11,'Text','','WordWrap',true);
            temp.grid   = uigridlayout(uipanel(f,'Position',[10 40 280 100],'BorderType','none'),'ColumnWidth',{'1x','1x'},'Scrollable','on');
            temp.labels = [];
            temp.edits  = [];
            drop        = uidropdown(f,'Position',[10 150 140 20],'BackgroundColor','w','Items',obj.blockList(:,1),'ItemsData',1:size(obj.blockList,1),'Value',find(ismember(obj.blockList(:,1),'Block_Population')),'Interruptible','off','BusyAction','cancel','UserData',temp,'ValueChangedFcn',@obj.cb_selectblock);
            uibutton(f,'Position',[140 8 70 24],'Text','Add','Interruptible','off','BusyAction','cancel','ButtonpushedFcn',{@obj.cb_createblock,drop});
            uibutton(f,'Position',[220 8 70 24],'Text','Cancel','Interruptible','off','BusyAction','cancel','ButtonpushedFcn',@(h,~)delete(h.Parent));
            obj.cb_selectblock(drop);
        end
        %% Select a block in the drop-down
        function cb_selectblock(obj,h,~)
            h.UserData.tip.Text = obj.blockList{h.Value,2};
            delete(h.UserData.labels);
            delete(h.UserData.edits);
            h.UserData.labels = [];
            h.UserData.edits  = [];
            Parameter = obj.blockList{h.Value,3};
            h.UserData.grid.RowHeight = repmat({22},1,max(1,size(Parameter,1)));
            if isempty(Parameter)
                h.UserData.labels = GUI.APP(1,[1 2],uilabel(h.UserData.grid,'Text','This block does not contain any parameter.','HorizontalAlignment','center'));
            else
                for i = 1 : size(Parameter,1)
                    h.UserData.labels = [h.UserData.labels,GUI.APP(i,1,uilabel(h.UserData.grid,'Text',Parameter{i,1},'Tooltip',Parameter{i,3}))];
                    h.UserData.edits  = [h.UserData.edits,GUI.APP(i,2,uieditfield(h.UserData.grid,'Value',Parameter{i,2},'Tooltip',Parameter{i,3},'HorizontalAlignment','right'))];
                end
            end
        end
        %% Create a new block
        function cb_createblock(obj,h,~,drop)
            index = drop.Value;
            para  = cell(1,size(obj.blockList{index,3},1));
            for i = 1 : length(para)
                para{i} = str2num(drop.UserData.edits(i).Value);
            end
            block     = feval(obj.blockList{index,1},para{:});
            app       = CreateBlockObject(obj,true);
            obj.Graph = addnode(obj.Graph,table(block,app));
            RefreshText(obj.Graph,obj.app.canvas.UserData.features,numnodes(obj.Graph),1);
            delete(h.Parent);
        end
        %% Move on the figure
        function cb_figmove(obj,~,~)
            h = obj.app.canvas.UserData.object;
            if isempty(h.UserData)
                % Move a block
                h.Children(2).Position(1:2) = h.Children(2).Position(1:2) + obj.app.canvas.CurrentPoint(1,1:2) - obj.app.canvas.UserData.lastpos;
                h.Children(1).Position(1:2) = h.Children(2).Position(1:2) + [7 4.5];
                obj.app.canvas.UserData.connectable = false;
                blockNo = find(h==obj.Graph.Nodes.app);
                % Move the lines starting from the block
                for i = outedges(obj.Graph,blockNo)'
                    obj.Graph.Edges.app(i).UserData = obj.Graph.Edges.app(i).UserData + 0.5*(obj.app.canvas.CurrentPoint(1,1:2)-obj.app.canvas.UserData.lastpos);
                    FlushCurve(obj.Graph.Edges.app(i),h.Children(1).Position(1:2),[obj.Graph.Edges.app(i).Children(2).XData(end),obj.Graph.Edges.app(i).Children(2).YData(end)]);
                end
                % Move the lines ending at the block
                for i = inedges(obj.Graph,blockNo)'
                    obj.Graph.Edges.app(i).UserData = obj.Graph.Edges.app(i).UserData + 0.5*(obj.app.canvas.CurrentPoint(1,1:2)-obj.app.canvas.UserData.lastpos);
                    FlushCurve(obj.Graph.Edges.app(i),[obj.Graph.Edges.app(i).Children(2).XData(1),obj.Graph.Edges.app(i).Children(2).YData(1)],h.Children(1).Position(1:2));
                end
            else
                % Move a line
                h.UserData = h.UserData + 2*(obj.app.canvas.CurrentPoint(1,1:2)-obj.app.canvas.UserData.lastpos);
                FlushCurve(h,[h.Children(2).XData(1),h.Children(2).YData(1)],[h.Children(2).XData(end),h.Children(2).YData(end)]);
            end
            obj.app.canvas.UserData.lastpos = obj.app.canvas.CurrentPoint(1,1:2);
        end
        %% Click on a block
        function cb_block(obj,h,event)
            if event.Button ~= 3
                % Left-click on the line
                if obj.app.canvas.UserData.connectable && obj.app.canvas.UserData.object ~= h
                    originNo = find(obj.Graph.Nodes.app==obj.app.canvas.UserData.object);
                    termiNo  = find(obj.Graph.Nodes.app==h);
                    if findedge(obj.Graph,originNo,termiNo) == 0
                        % Create a new line
                        obj.cb_checkdetail([originNo,termiNo]);
                        return;
                    end
                end
                rotate3d(obj.GUI.app.figure,'off');
                pan(obj.GUI.app.figure,'off');
                zoom(obj.GUI.app.figure,'off');
                set(obj.GUI.app.figure,'WindowButtonMotionFcn',@obj.cb_figmove,'WindowButtonUpFcn',@(h,~)set(h,'WindowButtonMotionFcn',@(~,~)[],'WindowButtonUpFcn',[]));
            else
                % Right-click on the block
                set(obj.app.menu.items(3:4),'Enable',true);
                obj.app.menu.show();
            end
            % Highlight the block
            HighlightObject(obj.app.canvas,h);
        end
        %% Click on a line
        function cb_line(obj,h,event)
            % Highlight the line
            HighlightObject(obj.app.canvas,h);
            if event.Button ~= 3
                % Left-click on the line
                rotate3d(obj.GUI.app.figure,'off');
                pan(obj.GUI.app.figure,'off');
                zoom(obj.GUI.app.figure,'off');
                set(obj.GUI.app.figure,'WindowButtonMotionFcn',@obj.cb_figmove,'WindowButtonUpFcn',@(h,~)set(h,'WindowButtonMotionFcn',@(~,~)[],'WindowButtonUpFcn',[]));
            else
                % Right-click on the line
                set(obj.app.menu.items(3:4),'Enable',false);
                obj.app.menu.show();
            end
        end
        %% Select the features to display
        function cb_selectfeature(obj,~,~)
            P    = obj.GUI.app.figure.Position;
            C    = obj.GUI.app.figure.CurrentPoint;
            f    = uifigure('Name','Select features to display','Position',[P(1)+C(1) P(2)+C(2)-200 300 180],'icon',obj.GUI.icon.logo1,'WindowStyle','modal','Resize','off');
            grid = uigridlayout(uipanel(f,'Position',[10 40 280 140],'BorderType','none'),'RowHeight',{20,20,20,20,20,20},'ColumnWidth',{'1x','1x'},'RowSpacing',2);
            GUI.APP(1,1,uilabel(grid,'Text','Feature of blocks','FontWeight','bold'));
            check(1) = GUI.APP(2,1,uicheckbox(grid,'Text','Block number','Value',obj.app.canvas.UserData.features(1)));
            check(2) = GUI.APP(2,2,uicheckbox(grid,'Text','Block name','Value',obj.app.canvas.UserData.features(2)));
            check(3) = GUI.APP(3,1,uicheckbox(grid,'Text','#Input solutions','Value',obj.app.canvas.UserData.features(3)));
            check(4) = GUI.APP(3,2,uicheckbox(grid,'Text','#Output solutions','Value',obj.app.canvas.UserData.features(4)));
            check(5) = GUI.APP(4,1,uicheckbox(grid,'Text','#Parameters','Value',obj.app.canvas.UserData.features(5)));
            check(6) = GUI.APP(4,2,uicheckbox(grid,'Text','#Training times','Value',obj.app.canvas.UserData.features(6)));
            GUI.APP(5,1,uilabel(grid,'Text','Feature of lines','FontWeight','bold'));
            check(7) = GUI.APP(6,1,uicheckbox(grid,'Text','Ratio of solutions','Value',obj.app.canvas.UserData.features(7)));
            check(8) = GUI.APP(6,2,uicheckbox(grid,'Text','Number of solutions','Value',obj.app.canvas.UserData.features(8)));
            uibutton(f,'Position',[140 8 70 24],'Text','OK','Interruptible','off','BusyAction','cancel','ButtonpushedFcn',{@obj.cb_showfeature,check});
            uibutton(f,'Position',[220 8 70 24],'Text','Cancel','Interruptible','off','BusyAction','cancel','ButtonpushedFcn',@(h,~)delete(h.Parent));
        end
        %% Change the features of all blocks and lines
        function cb_showfeature(obj,h,~,check)
            obj.app.canvas.UserData.features = [check.Value];
            RefreshText(obj.Graph,obj.app.canvas.UserData.features);
            delete(h.Parent);
        end
        %% Arrange the blocks
        function cb_arrange(obj,~,~)
            tempG  = obj.Graph;
            tempG  = rmedge(tempG,find(ismember(tempG.Edges.EndNodes,tempG.Edges.EndNodes(:,[2,1]),'rows')));
            h      = plot(obj.app.canvas,tempG,'Visible',false,'Layout','layered','Direction','right');
            xcoord = (h.XData-min(h.XData)+1e-6)./(max(h.XData)-min(h.XData)+2e-6)*80 + 3;
            ycoord = (h.YData-min(h.YData)+1e-6)./(max(h.YData)-min(h.YData)+2e-6)*80 + 6;
            delete(h);
            axis(obj.app.canvas,[0 100 0 100]);
            % Move the blocks
            for i = 1 : numnodes(obj.Graph)
                obj.Graph.Nodes.app(i).Children(2).Position(1:2) = [xcoord(i),ycoord(i)];
                obj.Graph.Nodes.app(i).Children(1).Position(1:2) = obj.Graph.Nodes.app(i).Children(2).Position(1:2) + [7 4.5];
            end
            % Move the lines
            for i = 1 : numedges(obj.Graph)
                [originNo,termiNo] = findedge(obj.Graph,i);
                obj.Graph.Edges.app(i).UserData = (obj.Graph.Nodes.app(originNo).Children(1).Position(1:2)+obj.Graph.Nodes.app(termiNo).Children(1).Position(1:2))/2;
                if findedge(obj.Graph,termiNo,originNo)
                    % Translate the centers of bidirectional edges
                    obj.Graph.Edges.app(i).UserData(2) = obj.Graph.Edges.app(i).UserData(2) + 16*(-1)^(originNo>termiNo);
                end
                FlushCurve(obj.Graph.Edges.app(i),obj.Graph.Nodes.app(originNo).Children(1).Position(1:2),obj.Graph.Nodes.app(termiNo).Children(1).Position(1:2));
            end
        end
        %% Show the details of a block or line
        function cb_checkdetail(obj,ui,~)
            if nargin > 2
                h = obj.app.canvas.UserData.object;
                if isempty(h)
                    return;
                elseif isempty(h.UserData)
                    blockNo = find(h==obj.Graph.Nodes.app);
                else
                    blockNo = obj.Graph.Edges.EndNodes(h==obj.Graph.Edges.app,:);
                end
            else
                h = [];
                blockNo = ui;
            end
            P = obj.GUI.app.figure.Position;
            C = obj.GUI.app.figure.CurrentPoint;
            if ~isempty(h) && isempty(h.UserData)
                % Show the details of a block
                f         = uifigure('Name','Details of block','Position',[P(1)+C(1) P(2)+C(2)-260 400 240],'icon',obj.GUI.icon.logo1,'WindowStyle','modal','Resize','off');
                index     = find(ismember(obj.blockList(:,1),class(obj.Graph.Nodes.block(blockNo))));
                Parameter = obj.blockList{index,3};
                uilabel(f,'Position',[15 210 370 30],'HorizontalAlignment','center','Text',[GetBlockName(blockNo,obj.Graph.Nodes.block(blockNo),[1 1]),' (',obj.blockList{index,2},')'],'BackgroundColor',h.Children(2).FaceColor,'FontWeight','bold','WordWrap',true);
                grid      = uigridlayout(uipanel(f,'Position',[10 40 380 170],'BorderType','none'),'RowHeight',[{28,28,18,18},repmat({28},1,size(Parameter,1))],'ColumnWidth',{'5x','5x'},'RowSpacing',2,'Scrollable','on');
                GUI.APP(1,1,uilabel(grid,'Text','Predecessors:'));
                GUI.APP(1,2,uilabel(grid,'Text',[strjoin(arrayfun(@(s)GetBlockName(s,obj.Graph.Nodes.block(s),[1 1]),predecessors(obj.Graph,blockNo),'UniformOutput',false),', '),' (',GetInOutNumber(obj.Graph,blockNo,1),' solutions)'],'WordWrap',true));
                GUI.APP(2,1,uilabel(grid,'Text','Successors:'));
                GUI.APP(2,2,uilabel(grid,'Text',[strjoin(arrayfun(@(s)GetBlockName(s,obj.Graph.Nodes.block(s),[1 1]),successors(obj.Graph,blockNo),'UniformOutput',false),', '),' (',GetInOutNumber(obj.Graph,blockNo,2),' solutions)'],'WordWrap',true));
                GUI.APP(3,1,uilabel(grid,'Text','Number of parameters:'));
                GUI.APP(3,2,uilabel(grid,'Text',num2str(length(obj.Graph.Nodes.block(blockNo).parameter))));
                GUI.APP(4,1,uilabel(grid,'Text','Number of training times:'));
                GUI.APP(4,2,uilabel(grid,'Text',num2str(obj.Graph.Nodes.block(blockNo).trainTime)));
                for i = 1 : size(Parameter,1)
                    GUI.APP(4+i,1,uilabel(grid,'Text',['<',Parameter{i,1},'> ',Parameter{i,3}],'WordWrap',true));
                    GUI.APP(4+i,2,uilabel(grid,'Text',num2str(obj.Graph.Nodes.block(blockNo).(Parameter{i,1}))));
                end
                uibutton(f,'Position',[320 8 70 24],'Text','OK','Interruptible','off','BusyAction','cancel','ButtonpushedFcn',@(h,~)delete(h.Parent));
            else
                % Show the details of a line
                f    = uifigure('Name','Details of line','Position',[P(1)+C(1) P(2)+C(2)-150 260 130],'icon',obj.GUI.icon.logo1,'WindowStyle','modal','Resize','off');
                grid = uigridlayout(uipanel(f,'Position',[10 40 240 90],'BorderType','none'),'RowHeight',{22,22,22},'ColumnWidth',{'5x','5x'},'RowSpacing',2,'Scrollable','on');
                GUI.APP(1,1,uilabel(grid,'Text','Origin'));
                GUI.APP(1,2,uilabel(grid,'Text',GetBlockName(blockNo(1),obj.Graph.Nodes.block(blockNo(1)),[1 1]),'HorizontalAlignment','right'));
                GUI.APP(2,1,uilabel(grid,'Text','Terminus'));
                GUI.APP(2,2,uilabel(grid,'Text',GetBlockName(blockNo(2),obj.Graph.Nodes.block(blockNo(2)),[1 1]),'HorizontalAlignment','right'));
                GUI.APP(3,1,uilabel(grid,'Text','Ratio of solutions'));
                if isempty(h)
                    edit = GUI.APP(3,2,uieditfield(grid,'numeric','Value',1,'Limits',[0 1]));
                else
                    edit = GUI.APP(3,2,uieditfield(grid,'numeric','Value',obj.Graph.Edges.Weight(h==obj.Graph.Edges.app),'Limits',[0 1]));
                end
                uibutton(f,'Position',[100 8 70 24],'Text','OK','Interruptible','off','BusyAction','cancel','ButtonpushedFcn',{@obj.cb_createline,blockNo,edit});
                uibutton(f,'Position',[180 8 70 24],'Text','Cancel','Interruptible','off','BusyAction','cancel','ButtonpushedFcn',@(h,~)delete(h.Parent));
            end
        end
        %% Create a new line
        function cb_createline(obj,h,~,blockNo,edit)
            lineNo = findedge(obj.Graph,blockNo(1),blockNo(2));
            if lineNo == 0
                if edit.Value > 0
                    % Add a line
                    app       = CreateLineObject(obj,blockNo,true);
                    obj.Graph = addedge(obj.Graph,blockNo(1),blockNo(2),table(edit.Value,app,'VariableNames',{'Weight','app'}));
                    RefreshText(obj.Graph,obj.app.canvas.UserData.features,find(app==obj.Graph.Edges.app),2);
                end
            else
                if edit.Value > 0
                    obj.Graph.Edges.Weight(lineNo) = edit.Value;
                    RefreshText(obj.Graph,obj.app.canvas.UserData.features,lineNo,2);
                else
                    % Delete an existing line
                    delete(obj.Graph.Edges.app(lineNo));
                    obj.Graph = rmedge(obj.Graph,lineNo);
                    obj.app.canvas.UserData.object = [];
                end
            end
            delete(h.Parent);
        end
        %% Change the color of a block or line
        function cb_setcolor(obj,~,~)
            if ~isempty(obj.app.canvas.UserData.object)
                h = obj.app.canvas.UserData.object;
                if isempty(h.UserData)
                    color = uisetcolor(h.Children(2).FaceColor);
                    if length(color) == 3
                        h.Children(2).FaceColor = color;
                    end
                else
                    color = uisetcolor(h.Children(2).UserData);
                    if length(color) == 3
                        h.Children(2).UserData = color;
                        HighlightObject(obj.app.canvas);
                    end
                end
                figure(obj.GUI.app.figure);
            end
        end
        %% Open the file of a block
        function cb_openfile(obj,~,~)
            if ~isempty(obj.app.canvas.UserData.object)
                blockNo = find(obj.app.canvas.UserData.object==obj.Graph.Nodes.app);
                strP    = which(class(obj.Graph.Nodes.block(blockNo))); % Name of pcode
                web(['file://',strP(1:end-1),'m'],'-browser');
            end
        end
        %% Delete a block or line
        function cb_delete(obj,~,~)
            h = obj.app.canvas.UserData.object;
            if ~isempty(h)
                if isempty(h.UserData)
                    % Delete a block
                    blockNo = find(h==obj.Graph.Nodes.app);
                    delete(h);
                    if numedges(obj.Graph) > 0
                        delete(obj.Graph.Edges.app(outedges(obj.Graph,blockNo)));
                        delete(obj.Graph.Edges.app(inedges(obj.Graph,blockNo)));
                    end
                    obj.Graph = rmnode(obj.Graph,blockNo);
                    RefreshText(obj.Graph,obj.app.canvas.UserData.features);
                else
                    % Delete a line
                    lineNo = find(h==obj.Graph.Edges.app);
                    delete(h);
                    obj.Graph = rmedge(obj.Graph,lineNo);
                end
                obj.app.canvas.UserData.object      = [];
                obj.app.canvas.UserData.connectable = false;
            end
        end
        %% Load a block
        function cb_loadblock(obj,~,~)
            [Name,Path] = uigetfile({'*.mat','MAT file'});
            figure(obj.GUI.app.figure);
            if ischar(Name)
                try
                    data    = struct2cell(load(fullfile(Path,Name),'-mat'));
                    success = false;
                    for i = 1 : length(data)
                        if isa(data{i},'BLOCK')
                            for j = 1 : length(data{i})
                                % Add a block object
                                block     = data{i}(j);
                                app       = CreateBlockObject(obj,true);
                                obj.Graph = addnode(obj.Graph,table(block,app));
                                RefreshText(obj.Graph,obj.app.canvas.UserData.features,numnodes(obj.Graph),1);
                                success   = true;
                            end
                        end
                    end
                    assert(success);
                catch
                    uialert(obj.GUI.app.figure,'The selected file does not contain BLOCK objects.','Error');
                end
            end
        end
        %% Load an algorithm
        function cb_loadalgorithm(obj,~,~)
            if numnodes(obj.Graph)==0 || strcmp('OK',uiconfirm(obj.GUI.app.figure,'Load a new algorithm? All existing blocks will be deleted.','Warning','Icon','warning'))
                [Name,Path] = uigetfile({'*.mat','MAT file'});
                figure(obj.GUI.app.figure);
                if ischar(Name)
                    try
                        load(fullfile(Path,Name),'Blocks','Graph','-mat');
                        SetNewAlgorithm(obj,Blocks,Graph);
                    catch
                        uialert(obj.GUI.app.figure,['Fail to load the variables Blocks and Graph from ',fullfile(Path,Name)],'Error');
                    end
                end
            end
        end
        %% Save a block
        function cb_saveblock(obj,~,~)
            if ~isempty(obj.app.canvas.UserData.object)
                [Name,Path] = uiputfile({'*.mat','MAT file'},'','myBlock');
                figure(obj.GUI.app.figure);
                if ischar(Name)
                    try
                        block = obj.Graph.Nodes.block(obj.app.canvas.UserData.object==obj.Graph.Nodes.app);
                        save(fullfile(Path,Name),'block','-mat');
                    catch err
                        uialert(obj.GUI.app.figure,'Fail to save the block, please refer to the command window for details.','Error');
                        rethrow(err);
                    end
                end
            end
        end
        %% Save an algorithm
        function cb_savealgorithm(obj,~,~,Name,Path)
            if nargin < 5
                [Name,Path] = uiputfile({'*.mat','MAT file'},'','myAlgorithm');
                figure(obj.GUI.app.figure);
            end
            if ischar(Name)
                try
                    Blocks = obj.Graph.Nodes.block';
                    Graph  = full(adjacency(obj.Graph,'weighted'));
                    save(fullfile(Path,Name),'Blocks','Graph','-mat');
                catch err
                    uialert(obj.GUI.app.figure,'Fail to save the algorithm, please refer to the command window for details.','Error');
                    rethrow(err);
                end
            end
        end
        %% Generate the code of the algorithm
        function cb_generatecode(obj,~,~)
            if numnodes(obj.Graph) > 0
                [Name,Path] = uiputfile({'*.m','MATLAB function'},'',fullfile(cd,'myAlgorithm'));
                figure(obj.GUI.app.figure);
                if ischar(Name)
                    % Save the source code
                    try
                        [~,name] = fileparts(Name);
                        Code = {['classdef ',name,' < GEA % < ALGORITHM']
                               ''
                               '    methods'
                               ['        function obj = ',name,'(varargin)']
                               '            obj    = obj@GEA(varargin{:});'};
                        for i = 1 : numnodes(obj.Graph)
                            Parameter = obj.blockList{ismember(obj.blockList(:,1),class(obj.Graph.Nodes.block(i))),3};
                            paras     = cell(1,size(Parameter,1));
                            for j = 1 : size(Parameter,1)
                                paras{j} = num2str(obj.Graph.Nodes.block(i).(Parameter{j,1}));
                            end
                            if i == 1
                                str = '            Blocks = [';
                            else
                                str = '                      ';
                            end
                            str = [str,class(obj.Graph.Nodes.block(i)),'(',strjoin(paras,','),')'];
                            if i < numnodes(obj.Graph)
                                str = [str,',...'];
                            else
                                str = [str,'];'];
                            end
                            Code = [Code;str];
                        end
                        adjMat = num2str(full(adjacency(obj.Graph,'weighted')));
                        for i = 1 : size(adjMat,1)
                            if i == 1
                                str = '            Graph = [';
                            else
                                str = '                     ';
                            end
                            str = [str,adjMat(i,:)];
                            if i >= size(adjMat,1)
                                str = [str,'];'];
                            end
                            Code = [Code;str];
                        end
                        Code = [Code;
                               {''}
                               ['            load(''',name,'.mat'',''Blocks'',''Graph'');']
                               {''}
                               '            obj.parameter = {Blocks,Graph};'
                               '        end'
                               '    end'
                               'end'];
                        fid  = fopen(fullfile(Path,Name),'wt');
                        for i = 1 : length(Code)
                            fprintf(fid,'%s\n',Code{i});
                        end
                        fclose(fid);
                        web(['file://',fullfile(Path,Name)],'-browser');
                    catch err
                        uialert(obj.GUI.app.figure,'Fail to generate code, please refer to the command window for details.','Error');
                        rethrow(err);
                    end
                    % Save the algorithm
                    obj.cb_savealgorithm([],[],[name,'.mat'],Path);
                end
            end
        end
        %% Valid the algorithm
        function cb_validation(obj,~,~)
            try
                % Generate the ALGORITHM object
                Blocks = obj.Graph.Nodes.block;
                Graph  = adjacency(obj.Graph,'weighted');
                Blocks.Validity(Graph);
                ALG = GEA('parameter',{Blocks,Graph},'outputFcn',@(~,~)[]);
                % Generate the PROBLEM object
                [name,para] = GUI.GetParameterSetting(obj.app.listD.items(1));
                PRO = feval(name,'N',para{1},'M',para{2},'D',para{3},'maxFE',para{1}+1,'parameter',para(5:end));
                % Execute the algorithm for one generation
                ALG.Solve(PRO);
                RefreshText(obj.Graph,obj.app.canvas.UserData.features);
                uialert(obj.GUI.app.figure,'The algorithm is valid.','Valid algorithm','Icon','success');
            catch err
                err = addCause(err,MException('','The algorithm is invalid'));
                uialert(obj.GUI.app.figure,sprintf('%s, since %s',err.cause{end}.message,err.message),'Invalid algorithm');
                rethrow(err);
            end
        end
        %% Select a predefined algorithm
        function cb_selectAlgorithm(obj,h,~)
            if h.Value == 1
                if ~isempty(h.UserData)
                    SetNewAlgorithm(obj);
                    obj.Graph = h.UserData;
                    if numnodes(obj.Graph) > 0
                        set(obj.Graph.Nodes.app,'Visible',true);
                    end
                    if numedges(obj.Graph) > 0
                        set(obj.Graph.Edges.app,'Visible',true);
                    end
                    RefreshText(obj.Graph,obj.app.canvas.UserData.features);
                    h.UserData = [];
                end
            else
                switch h.Value
                    case 2
                        Blocks = [Block_Population(),Block_Tournament(200,3),Block_Crossover(2,3),Block_Mutation(3),Block_Selection(100)];
                        Graph  = [0 1 0 0 1
                                  0 0 1 0 0
                                  0 0 0 1 0
                                  0 0 0 0 1
                                  1 0 0 0 0];
                    case 3
                        Blocks = [Block_Population(),Block_Tournament(100,10),Block_Crossover(2,3),Block_Mutation(3),Block_Selection(100),...
                                  Block_Population(),Block_Tournament(100,10),Block_Crossover(2,3),Block_Mutation(3),Block_Selection(100)];
                        Graph  = [0 1 0 0 1 0 0 0 0 0
                                  0 0 1 0 0 0 0 0 0 0
                                  0 0 0 1 0 0 0 0 0 0
                                  0 0 0 0 1 0 0 0 0 1
                                  1 0 0 0 0 0 0 0 0 0
                                  0 0 0 0 0 0 1 0 0 1
                                  0 0 0 0 0 0 0 1 0 0
                                  0 0 0 0 0 0 0 0 1 0
                                  0 0 0 0 1 0 0 0 0 1
                                  0 0 0 0 0 1 0 0 0 0];
                    case 4
                        Blocks = [Block_Population(),Block_Tournament(200,10),Block_Tournament(200,10),Block_Tournament(200,10),...
                                  Block_Exchange(3),Block_Exchange(3),Block_Exchange(3),Block_Exchange(3),...
                                  Block_Crossover(2,3),Block_Mutation(3),Block_Selection(100)];
                        Graph  = [0 1 1 1 0 0 0 0 0 0 1
                                  0 0 0 0 1/4 1/4 1/4 1/4 0 0 0
                                  0 0 0 0 1/4 1/4 1/4 1/4 0 0 0
                                  0 0 0 0 1/4 1/4 1/4 1/4 0 0 0
                                  0 0 0 0 0 0 0 0 1 0 0
                                  0 0 0 0 0 0 0 0 1 0 0
                                  0 0 0 0 0 0 0 0 1 0 0
                                  0 0 0 0 0 0 0 0 1 0 0
                                  0 0 0 0 0 0 0 0 0 1 0
                                  0 0 0 0 0 0 0 0 0 0 1
                                  1 0 0 0 0 0 0 0 0 0 0];
                end
                if isempty(h.UserData)
                    if numnodes(obj.Graph) > 0
                        set(obj.Graph.Nodes.app,'Visible',false);
                    end
                    if numedges(obj.Graph) > 0
                        set(obj.Graph.Edges.app,'Visible',false);
                    end
                    h.UserData = obj.Graph;
                    obj.Graph  = digraph();
                end
                SetNewAlgorithm(obj,Blocks,Graph);
            end
        end
        %% Update the problems in the list
        function cb_filter(obj,~,~,index)
            GUI.UpdateAlgProList(index,obj.app.stateC,obj.app.listC,obj.app.dropC,obj.app.labelC(5),obj.GUI.proList);
        end
        %% Determine the file to load and save training results
        function cb_filepath(obj,~,~,edit)
            [Name,Path] = uigetfile({'*.mat','MAT file'},'',fileparts(edit.Value));
            figure(obj.GUI.app.figure);
            if ischar(Name)
                edit.Value = fullfile(Path,Name);
            end
        end
        %% Start the testing
        function cb_test(obj,~,~)
            if strcmp(obj.app.buttonE(1).Text,'Pause')
                obj.app.buttonE(1).Text = 'Continue';
            elseif strcmp(obj.app.buttonE(1).Text,'Continue')
                obj.app.buttonE(1).Text = 'Pause';
            else
                % Generate the ALGORITHM object
                try
                    Blocks = obj.Graph.Nodes.block;
                    Graph  = adjacency(obj.Graph,'weighted');
                    Blocks.Validity(Graph);
                    ALG = GEA('parameter',{Blocks,Graph},'outputFcn',@obj.outputFcnTest,'save',1);
                catch err
                    err = addCause(err,MException('','The algorithm is invalid'));
                    uialert(obj.GUI.app.figure,sprintf('%s, since %s',err.cause{end}.message,err.message),'Invalid algorithm');
                    rethrow(err);
                end
                % Generate the PROBLEM object
                try
                    [name,para] = GUI.GetParameterSetting(obj.app.listD.items(1));
                    PRO = feval(name,'N',para{1},'M',para{2},'D',para{3},obj.app.listD.items(1).label(4).Text,para{4},'parameter',para(5:end));
                catch err
                    uialert(obj.GUI.app.figure,err.message,'Invalid parameter settings');
                    return;
                end
                % Update the data
                obj.data = {ALG,PRO};
                % Update the GUI
                set([obj.GUI.app.button,obj.app.buttonA,obj.app.dropA,obj.app.stateC,obj.app.listC,obj.app.dropC,obj.app.editD,obj.app.checkD,obj.app.buttonD],'Enable',false);
                obj.app.listD.Enable    = false;
                obj.app.buttonE(1).Text = 'Pause';
                set(obj.app.buttonE(2:3),'Enable',true);
                delete(obj.app.axesE.Children);
                % Execute the algorithm
                try
                    ALG.Solve(PRO);
                    obj.cb_stoptest();
                catch err
                    uialert(obj.GUI.app.figure,'The algorithm terminates unexpectedly, please refer to the command window for details.','Error');
                    obj.cb_stoptest();
                    rethrow(err);
                end
            end
        end
        %% Stop the testing
        function cb_stoptest(obj,~,~,metValue)
            if nargin < 4
                RefreshText(obj.Graph,obj.app.canvas.UserData.features);
                set([obj.GUI.app.button,obj.app.buttonA,obj.app.dropA,obj.app.stateC,obj.app.listC,obj.app.dropC,obj.app.editD,obj.app.buttonD(1:2)],'Enable',true);
                obj.app.checkD.Enable     = ~isempty(ver('parallel'));
                obj.app.listD.Enable      = true;
                obj.app.buttonE(1).Text   = 'Start';
                obj.app.buttonE(2).Enable = false;
                obj.app.buttonE(3).Enable = ~isempty(obj.data{1}.result);
            end
            if isempty(obj.data{1}.result)
                obj.data = {};
            else
                if nargin < 4
                    if obj.data{2}.M == 1
                        metValue = obj.data{2}.CalMetric('Min_value',obj.data{1}.result{end});
                    else
                        metValue = obj.data{2}.CalMetric('HV',obj.data{1}.result{end});
                    end
                end
                Draw(obj.app.axesE);
                if obj.data{2}.M == 1
                    obj.data{2}.DrawDec(obj.data{1}.result{end});
                else
                    obj.data{2}.DrawObj(obj.data{1}.result{end});
                end
                obj.app.axesE.FontSize = 11;
                if obj.data{2}.M == 1
                    obj.app.labelE.Text = sprintf('Min value: %.2e',metValue);
                else
                    obj.app.labelE.Text = sprintf('HV: %.2e',metValue);
                end
            end
        end
        %% Start the training
        function cb_train(obj,~,~)
            if strcmp(obj.app.buttonD(2).Text,'Pause')
                obj.app.buttonD(2).Text = 'Continue';
            elseif strcmp(obj.app.buttonD(2).Text,'Continue')
                obj.app.buttonD(2).Text = 'Pause';
            else
                % Generate the ALGORITHM object of algorithm
                try
                    Blocks = obj.Graph.Nodes.block;
                    Graph  = adjacency(obj.Graph,'weighted');
                    Blocks.Validity(Graph);
                    obj.data{1} = GEA('parameter',{Blocks,Graph},'outputFcn',@(~,~)[],'save',1);
                catch err
                    err = addCause(err,MException('','The algorithm is invalid'));
                    uialert(obj.GUI.app.figure,sprintf('%s, since %s',err.cause{end}.message,err.message),'Invalid algorithm');
                    rethrow(err);
                end
                % Generate the PROBLEM object of algorithm
                try
                    [name,para] = GUI.GetParameterSetting(obj.app.listD.items(1));
                    obj.data{2} = feval(name,'N',para{1},'M',para{2},'D',para{3},obj.app.listD.items(1).label(4).Text,para{4},'parameter',para(5:end));
                catch err
                    uialert(obj.GUI.app.figure,err.message,'Invalid parameter settings');
                    return;
                end
                % Update the GUI
                set([obj.GUI.app.button,obj.app.buttonA,obj.app.dropA,obj.app.stateC,obj.app.listC,obj.app.dropC,obj.app.editD,obj.app.checkD,obj.app.buttonD(1),obj.app.buttonE],'Enable',false);
                obj.app.listD.Enable      = false;
                obj.app.buttonD(2).Text   = 'Pause';
                obj.app.buttonD(3).Enable = true;
                delete(obj.app.axesD.Children);
                delete(obj.app.axesE.Children);
                try
                    % Generate the ALGORITHM object of trainer
                    obj.dataTrain{1} = GA('outputFcn',@obj.outputFcnTrain,'save',1);
                    % Generate the PROBLEM object of trainer
                    obj.dataTrain{2} = UserProblem('N',obj.app.editD(1).Value,'maxFE',obj.app.editD(2).Value,'D',length(Blocks.parameters),'lower',Blocks.lowers,'upper',Blocks.uppers,'initFcn',@obj.trainInit,'objFcn',@obj.trainObj,'once',true);
                    % Execute the training
                    obj.dataTrain{1}.Solve(obj.dataTrain{2});
                    obj.cb_stoptrain();
                catch err
                    obj.cb_stoptrain();
                    if ~strcmp(err.identifier,'PlatEMO:Termination')
                        uialert(obj.GUI.app.figure,'The trainer terminates unexpectedly, please refer to the command window for details.','Error');
                        rethrow(err);
                    end
                end
            end
        end
        %% Stop the training
        function cb_stoptrain(obj,~,~)
            RefreshText(obj.Graph,obj.app.canvas.UserData.features);
            set([obj.GUI.app.button,obj.app.buttonA,obj.app.dropA,obj.app.stateC,obj.app.listC,obj.app.dropC,obj.app.editD,obj.app.buttonD(1),obj.app.buttonE(1:2)],'Enable',true);
            obj.app.checkD.Enable     = ~isempty(ver('parallel'));
            obj.app.listD.Enable      = true;
            obj.app.buttonD(2).Text   = 'Start';
            obj.app.buttonD(3).Enable = false;
            obj.app.buttonE(3).Enable = ~isempty(obj.data{1}.result);
            if isempty(obj.data{1}.result)
                obj.data = {};
            end
            if isempty(obj.dataTrain{1}.result)
                obj.dataTrain = {};
            else
                if obj.data{2}.M == 1
                    obj.app.labelD.Text = sprintf('Min value: %.2e',min(obj.dataTrain{1}.result{end}.objs));
                else
                    obj.app.labelD.Text = sprintf('HV: %.2e',-min(obj.dataTrain{1}.result{end}.objs));
                end
            end
        end
        %% Output function of single test
        function outputFcnTest(obj,Algorithm,Problem)
            assert(strcmp(obj.app.buttonE(2).Enable,'on'),'PlatEMO:Termination','');
            obj.app.labelE.Text = sprintf('%.2f%%',Problem.FE/Problem.maxFE*100);
            if strcmp(obj.app.buttonE(1).Text,'Continue')
                waitfor(obj.app.buttonE(1),'Text');
            end
            assert(strcmp(obj.app.buttonE(2).Enable,'on'),'PlatEMO:Termination','');
        end
        %% Output function of training
        function outputFcnTrain(obj,Algorithm,Problem)
        	assert(strcmp(obj.app.buttonD(3).Enable,'on'),'PlatEMO:Termination','');
            % Show the results
            obj.app.labelD.Text = sprintf('%.2f%%',Problem.FE/Problem.maxFE*100);
            if obj.data{2}.M == 1
                point = [Algorithm.result{end,1},min(Algorithm.result{end,2}.objs)];
            else
                point = [Algorithm.result{end,1},-min(Algorithm.result{end,2}.objs)];
            end
            if isempty(obj.app.axesD.Children)
                plot(obj.app.axesD,point(1),point(2),'-k.','LineWidth',1);
                obj.app.axesD.YScale = 'linear';
            else
                [obj.app.axesD.Children(1).XData,obj.app.axesD.Children(1).YData] = deal([obj.app.axesD.Children(1).XData,point(1)],[obj.app.axesD.Children(1).YData,point(2)]);
            end
            if min(obj.app.axesD.Children(1).YData) > 0 && max(obj.app.axesD.Children(1).YData) > 10*min(obj.app.axesD.Children(1).YData)
                obj.app.axesD.YScale = 'log';
            end
            % Save the results
            try
                Population = Algorithm.result{end};
                [~,best]   = min(Population.objs);
                obj.Graph.Nodes.block.ParameterSet(Population(best).dec);
                Blocks = obj.Graph.Nodes.block';
                Graph  = full(adjacency(obj.Graph,'weighted'));
                save(obj.app.editD(4).Value,'Blocks','Graph','Population','-mat');
            catch err
                uialert(obj.GUI.app.figure,'Fail to save the population, please refer to the command window for details.','Error');
                rethrow(err);
            end
            if strcmp(obj.app.buttonD(2).Text,'Continue')
                waitfor(obj.app.buttonD(2),'Text');
            end
            assert(strcmp(obj.app.buttonD(3).Enable,'on'),'PlatEMO:Termination','');
        end
        %% Initialization function for training
        function PopDec = trainInit(obj,N)
            D = length(obj.Graph.Nodes.block.parameters);
            if exist(obj.app.editD(4).Value,'file')
                try
                    load(obj.app.editD(4).Value,'-mat','Population');
                    PopDec = Population.decs;
                catch
                    error(['Fail to load the variable Population from ',obj.app.editD(4).Value]);
                end
                assert(D==size(PopDec,2),'The current algorithm has %d parameters while the loaded Population has %d variables',D,size(PopDec,2));
            else
                PopDec = [];
            end
            if size(PopDec,1) > N
                PopDec = PopDec(1:N,:);
            elseif size(PopDec,1) < N
                PopDec = [PopDec;unifrnd(repmat(obj.Graph.Nodes.block.lowers,N-size(PopDec,1),1),repmat(obj.Graph.Nodes.block.uppers,N-size(PopDec,1),1))];
            end
        end
        %% Objective function for training
        function PopObj = trainObj(obj,PopDec)
            PopObj = zeros(size(PopDec,1),obj.app.editD(3).Value);
            if ~obj.app.checkD.Value    % Evaluation in sequence
                for i = 1 : size(PopObj,1)
                    obj.data{1}.parameter{1}.ParameterSet(PopDec(i,:));
                    for j = 1 : size(PopObj,2)
                        assert(strcmp(obj.app.buttonD(3).Enable,'on'),'PlatEMO:Termination','');
                        obj.data{1}.Solve(obj.data{2});
                        if obj.data{2}.M == 1
                            PopObj(i,j) = obj.data{2}.CalMetric('Min_value',obj.data{1}.result{end});
                            obj.cb_stoptest([],[],PopObj(i,j));
                        else
                            PopObj(i,j) = -obj.data{2}.CalMetric('HV',obj.data{1}.result{end});
                            obj.cb_stoptest([],[],-PopObj(i,j));
                        end
                        if strcmp(obj.app.buttonD(2).Text,'Continue')
                            waitfor(obj.app.buttonD(2),'Text');
                        end
                        assert(strcmp(obj.app.buttonD(3).Enable,'on'),'PlatEMO:Termination','');
                    end
                end
            else                        % Evaluation in parallel
                for b = obj.Graph.Nodes.block'
                    b.trainTime = b.trainTime + size(PopDec,1);
                end
                for i = 1 : numel(PopObj)
                    Future(i) = parfeval(@parallelFcn,1,obj.data{1},obj.data{2},PopDec(mod(i-1,size(PopDec,1))+1,:));
                end
                while ~all([Future.Read])
                    drawnow('limitrate');
                    if strcmp(obj.app.buttonD(3).Enable,'off')
                        cancel(Future);
                        error('PlatEMO:Termination','');
                    end
                    [r,ALG] = fetchNext(Future,0.01);
                    if ~isempty(r)
                        obj.data{1} = ALG;
                        [i,j]       = ind2sub(size(PopObj),r);
                        if obj.data{2}.M == 1
                            PopObj(i,j) = obj.data{2}.CalMetric('Min_value',obj.data{1}.result{end});
                            obj.cb_stoptest([],[],PopObj(i,j));
                        else
                            PopObj(i,j) = -obj.data{2}.CalMetric('HV',obj.data{1}.result{end});
                            obj.cb_stoptest([],[],-PopObj(i,j));
                        end
                    end
                    if strcmp(obj.app.buttonD(2).Text,'Continue')
                        waitfor(obj.app.buttonD(2),'Text');
                    end
                    if strcmp(obj.app.buttonD(3).Enable,'off')
                        cancel(Future);
                        error('PlatEMO:Termination','');
                    end
                end
            end
            PopObj = mean(PopObj,2) + std(PopObj,0,2);
        end
    end
end

%% Create the graphic objects of a block
function app = CreateBlockObject(obj,highlight)
    app = hggroup(obj.app.canvas,'ButtonDownFcn',@obj.cb_block);
    rectangle(app,'Position',[0 91 14 9],'FaceColor',[.95 .95 1],'EdgeColor',[.8 .8 .8],'Curvature',[0.1 0.1],'HitTest',false);
    text(app,7,95.5,'','FontSize',11,'HorizontalAlignment','center','Clipping',true,'HitTest',false);
    if highlight
        HighlightObject(obj.app.canvas,app);
    else
        HighlightObject(obj.app.canvas);
    end
end

%% Create the graphic objects of a line
function app = CreateLineObject(obj,blockNo,highlight)
    origin = obj.Graph.Nodes.app(blockNo(1)).Children(1).Position;
    termi  = obj.Graph.Nodes.app(blockNo(2)).Children(1).Position;
    app    = hggroup(obj.app.canvas,'ButtonDownFcn',@obj.cb_line,'UserData',(origin(1:2)+termi(1:2))/2);
    line(app,0,0,'Color','k','UserData',[0 0 0]);
    text(app,0,0,'>','Margin',0.1,'Color','k','BackgroundColor','w','FontSize',11,'HorizontalAlignment','center','Clipping',true,'HitTest',false);
    FlushCurve(app,origin(1:2),termi(1:2));
    obj.app.canvas.Children = obj.app.canvas.Children([2:end,1]);
    if highlight
        HighlightObject(obj.app.canvas,app);
    else
        HighlightObject(obj.app.canvas);
    end
end

%% Highlight a block or line in the canvas
function HighlightObject(canvas,app)
    if nargin > 1
        HighlightObject(canvas);
        app.Children(1).Color = 'b';
        if isempty(app.UserData)
            % Highlight a hggroup of block
            app.Children(2).EdgeColor   = 'b';
            canvas.UserData.connectable = true;
        else
            % Highlight a hggroup of line
            app.Children(2).Color = 'b';
        end
        canvas.UserData.object  = app;
        canvas.UserData.lastpos = canvas.CurrentPoint(1,1:2);
    elseif ~isempty(canvas.UserData.object)
        if isvalid(canvas.UserData.object)
            if isempty(canvas.UserData.object.UserData)
                % Unhighlight a hggroup of block
                canvas.UserData.object.Children(1).Color     = 'k';
                canvas.UserData.object.Children(2).EdgeColor = [.8 .8 .8];
            else
                % Unhighlight a hggroup of line
                set(canvas.UserData.object.Children,'Color',canvas.UserData.object.Children(2).UserData)
            end
        end
        canvas.UserData.object      = [];
        canvas.UserData.connectable = false;
    end
end

%% Refresh a Bezier curve
function FlushCurve(curve,origin,termi)
    mid = curve.UserData;
    t   = linspace(0,1,max(20,norm(termi-origin)/2))';
    P   = (1-t).^2.*origin + 2*t.*(1-t).*mid + t.^2.*termi;
    if P(1,1) < P(end,1) && curve.Children(1).String(1)=='<'
        curve.Children(1).String = [curve.Children(1).String(2:end),'>'];
    elseif P(1,1) >= P(end,1) && curve.Children(1).String(end)=='>'
        curve.Children(1).String = ['<',curve.Children(1).String(1:end-1)];
    end
    curve.Children(2).XData = P(:,1);
    curve.Children(2).YData = P(:,2);
    curve.Children(1).Rotation      = atan((P(ceil(end/2)+1,2)-P(ceil(end/2),2))./(P(ceil(end/2)+1,1)-P(ceil(end/2),1)-1e-6).*curve.Parent.Position(4)./curve.Parent.Position(3))/pi*180;
    curve.Children(1).Position(1:2) = P(ceil(end/2),:);
end

%% Refresh the text of a block or line
function RefreshText(Graph,features,no,type)
    if nargin < 3
        % Refresh the texts of all blocks and lines
        for i = 1 : numnodes(Graph)
            RefreshText(Graph,features,i,1);
        end
        for i = 1 : numedges(Graph)
            RefreshText(Graph,features,i,2);
        end
    elseif type == 1
        % Refresh the text of a block
        app   = Graph.Nodes.app(no);
        block = Graph.Nodes.block(no);
        str   = {};
        str   = [str,GetBlockName(no,block,features(1:2))];
        if features(3) && features(4)
            str = [str,['In: ',GetInOutNumber(Graph,no,1),' Out: ',GetInOutNumber(Graph,no,2)]];
        elseif features(3) && ~features(4)
            str = [str,['In: ',GetInOutNumber(Graph,no,1)]];
        elseif ~features(3) && features(4)
            str = [str,['Out: ',GetInOutNumber(Graph,no,2)]];
        end
        if features(5) && features(6)
            str = [str,['Para: ',num2str(length(block.parameter)),' Train: ',num2str(block.trainTime)]];
        elseif features(5) && ~features(6)
            str = [str,['Para: ',num2str(length(block.parameter))]];
        elseif ~features(5) && features(6)
            str = [str,['Train: ',num2str(block.trainTime)]];
        end
        app.Children(1).String = str;
    elseif type == 2
        % Refresh the text of a line
        app = Graph.Edges.app(no);
        if features(7) && features(8)
            str = [' ',num2str(Graph.Edges.Weight(no),3),'(',GetInOutNumber(Graph,no,3),') '];
        elseif features(7) && ~features(8)
            str = [' ',num2str(Graph.Edges.Weight(no),3),' '];
        elseif ~features(7) && features(8)
            str = [' (',GetInOutNumber(Graph,no,3),') '];
        else
            str = [];
        end
        if app.Children(2).XData(1) < app.Children(2).XData(end)
            str = [str,'>'];
        else
            str = ['<',str];
        end
        app.Children(1).String = str;
    end
end

%% Get the name of a block
function name = GetBlockName(no,block,type)
    str = class(block);
    if type(1) && type(2)
        name = ['#',num2str(no),' ',str(7:end)];
    elseif type(1) && ~type(2)
        name = ['#',num2str(no)];
    elseif ~type(1) && type(2)
        name = str(7:end);
    else
        name = [];
    end
end

%% Get the string of number of solutions of a block or line
function num = GetInOutNumber(Graph,no,type)
    function len = GetLen(out)
        if isa(out,'SOLUTION')
            len = length(out);
        else
            len = size(out,1);
        end
    end
    if type == 1
        % Get the number of input solutions of a block
        if no <= numnodes(Graph) 
            [eid,nid] = inedges(Graph,no);
            num = sum(arrayfun(@(e,n)floor(GetLen(Graph.Nodes.block(n).output)*Graph.Edges.Weight(e)),eid,nid));
        else
            num = 0;
        end
    elseif type == 2
        % Get the number of output solutions of a block
        num = GetLen(Graph.Nodes.block(no).output);
    elseif type == 3
        % Get the number of passed solutions of a line
        num = floor(GetLen(Graph.Nodes.block(Graph.Edges.EndNodes(no,1)).output)*Graph.Edges.Weight(no));
    end
    num = num2str(num);
end

%% Delete the current algorithm and set a new algorithm
function SetNewAlgorithm(obj,Blocks,Graph)
    if numnodes(obj.Graph) > 0
        delete(obj.Graph.Nodes.app);
    end
    if numedges(obj.Graph) > 0
        delete(obj.Graph.Edges.app);
    end
    if nargin > 1
        obj.Graph = digraph(Graph);
        obj.Graph.Nodes.block = Blocks(:);
        for i = 1 : numnodes(obj.Graph)
            obj.Graph.Nodes.app(i) = CreateBlockObject(obj,false);
        end
        for i = 1 : numedges(obj.Graph)
            obj.Graph.Edges.app(i) = CreateLineObject(obj,obj.Graph.Edges.EndNodes(i,:),false);
        end
        RefreshText(obj.Graph,obj.app.canvas.UserData.features);
        obj.cb_arrange();
    end
end

%% Function for parallelization of algorithm evaluation in training
function Algorithm = parallelFcn(Algorithm,Problem,Dec)
    Algorithm.parameter{1}.ParameterSet(Dec);
    Algorithm.Solve(Problem);
end