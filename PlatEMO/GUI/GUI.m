classdef GUI < handle
%GUI - The class of the main figure window of the platform.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        app = struct();     % All the components
        iconFolder;         % Folder of icon files
        algList;            % Algorithm list
        proList;            % Problem list
        metList;            % Metric list
    end
    methods
        %% Establish the figure window
        function obj = GUI()
            % Create the window
            obj.iconFolder   = fullfile(fileparts(mfilename('fullpath')),'icons');
            obj.app.figure   = uifigure('Name','PlatEMO v3.2','Position',[0 0 1200 650],'Interruptible','off','icon',fullfile(obj.iconFolder,'logo1.png'),'BusyAction','cancel','Visible','off','WindowButtonMotionFcn',@(~,~)[]);
            obj.app.maingrid = uigridlayout(obj.app.figure,'RowHeight',{25,80,'1x'},'ColumnWidth',{'1x'},'Padding',[0 0 0 0],'RowSpacing',0);
            
            % Create the tab buttons
            obj.app.grid(1)    = GUI.APP(1,1,uigridlayout(obj.app.maingrid,'RowHeight',{'1x'},'ColumnWidth',{80,80,'1x',100},'Padding',[0 0 0 0],'ColumnSpacing',0,'BackgroundColor',[0 .25 .45]));
            tempPanel          = GUI.APP(1,1,uipanel(obj.app.grid(1),'BorderType','none','BackgroundColor',[0 .25 .45]));
            obj.app.buttonT(1) = uibutton(tempPanel,'Position',[-5 -5 90 35],'Text','Modules','FontSize',14,'FontColor','w','BackgroundColor',[0 .25 .45],'ButtonpushedFcn',{@obj.cb_tab,1});
            tempPanel          = GUI.APP(1,2,uipanel(obj.app.grid(1),'BorderType','none','BackgroundColor',[0 .25 .45]));
            obj.app.buttonT(2) = uibutton(tempPanel,'Position',[-5 -5 90 35],'Text','Support','FontSize',14,'FontColor','w','BackgroundColor',[0 .25 .45],'ButtonpushedFcn',{@obj.cb_tab,2});
            tempPanel          = GUI.APP(1,4,uipanel(obj.app.grid(1),'BorderType','none','BackgroundColor',[0 .25 .45]));
            tempImage          = uiimage(tempPanel,'Position',[0 3 99 21],'ImageSource',fullfile(obj.iconFolder,'bar.png'));
            tempPanel2         = uipanel(tempPanel,'Position',[13 4 18 18],'BorderType','none','BackgroundColor',[.549 .6627 .7529]);
            obj.app.buttonT(3) = uibutton(tempPanel2,'Position',[-2 -2 24 24],'Text','','Icon',fullfile(obj.iconFolder,'researchgate.png'),'BackgroundColor',[.549 .6627 .7529],'Tooltip','ResearchGate','ButtonpushedFcn',@(~,~)web('https://www.researchgate.net/publication/312061859','-browser'));
            tempPanel2         = uipanel(tempPanel,'Position',[43 4 18 18],'BorderType','none','BackgroundColor',[.549 .6627 .7529]);
            obj.app.buttonT(4) = uibutton(tempPanel2,'Position',[-2 -2 24 24],'Text','','Icon',fullfile(obj.iconFolder,'github.png'),'BackgroundColor',[.549 .6627 .7529],'Tooltip','GitHub','ButtonpushedFcn',@(~,~)web('https://github.com/BIMK/PlatEMO','-browser'));
            tempPanel2         = uipanel(tempPanel,'Position',[73 4 18 18],'BorderType','none','BackgroundColor',[.549 .6627 .7529]);
            obj.app.buttonT(5) = uibutton(tempPanel2,'Position',[-2 -2 24 24],'Text','','Icon',fullfile(obj.iconFolder,'qq.png'),'BackgroundColor',[.549 .6627 .7529],'Tooltip','QQ','ButtonpushedFcn',@(~,~)web('https://jq.qq.com/?_wv=1027&k=KTo7jITo','-browser'));
            
            % Create the menu
            obj.app.grid(2)   = GUI.APP(2,1,uigridlayout(obj.app.maingrid,'RowHeight',{'1x',13,1},'ColumnWidth',{1,75,75,75,'1x',13,1},'Padding',[0 0 0 5],'RowSpacing',5));
            obj.app.button(1) = GUI.APP([1 2],2,uibutton(obj.app.grid(2),'Text',{'Test','Module'},'VerticalAlignment','bottom','FontSize',11,'Icon',fullfile(obj.iconFolder,'test.png'),'IconAlignment','top','Tooltip',{'Test one algorithm on a problem with specified parameter settings.','You can analyse the result and study the performance of the algorithm from various aspects.'},'ButtonpushedFcn',{@obj.cb_module,1}));
            obj.app.button(2) = GUI.APP([1 2],3,uibutton(obj.app.grid(2),'Text',{'Application','Module'},'VerticalAlignment','bottom','FontSize',11,'Icon',fullfile(obj.iconFolder,'application.png'),'IconAlignment','top','Tooltip',{'Use algorithms to solve your own problem.','You can design your own problem and solve it by the suggested algorithms.'},'ButtonpushedFcn',{@obj.cb_module,2}));
            obj.app.button(3) = GUI.APP([1 2],4,uibutton(obj.app.grid(2),'Text',{'Experiment','Module'},'VerticalAlignment','bottom','FontSize',11,'Icon',fullfile(obj.iconFolder,'experiment.png'),'IconAlignment','top','Tooltip',{'Do experiment on multiple algorithms and problems.','You can observe the statistical results shown in a table and save it as an Excel or LaTeX table.'},'ButtonpushedFcn',{@obj.cb_module,3}));
            obj.app.button(4) = GUI.APP([1 2],2,uibutton(obj.app.grid(2),'Visible',false,'Text',{'About','PlatEMO'},'VerticalAlignment','bottom','FontSize',11,'Icon',fullfile(obj.iconFolder,'author.png'),'IconAlignment','top','ButtonpushedFcn',@obj.cb_author));
            obj.app.button(5) = GUI.APP([1 2],3,uibutton(obj.app.grid(2),'Visible',false,'Text',{'User','Manual'},'VerticalAlignment','bottom','FontSize',11,'Icon',fullfile(obj.iconFolder,'help.png'),'IconAlignment','top','ButtonpushedFcn',@(~,~)web(['file://',fullfile(fileparts(fileparts(obj.iconFolder)),'manual.pdf')],'-browser')));
            obj.app.tip       = GUI.APP(2,6,uiimage(obj.app.grid(2),'ImageSource',fullfile(obj.iconFolder,'tip2.png'),'ImageClickedFcn',@obj.cb_fold));
            tempLine          = GUI.APP(3,[1 7],uipanel(obj.app.grid(2),'BackgroundColor',[.8 .8 .8]));
            
            % Create the modules
            movegui(obj.app.figure,'center');
            obj.app.figure.addlistener('CurrentPoint','PostSet',@obj.cb_motion);
            obj.readList();
            obj.cb_module([],[],1);
            obj.app.figure.Visible = 'on';
        end
    end
	methods(Access = private)
        %% Change the menu buttons
        function cb_tab(obj,~,~,type)
            switch type
                case 1
                    [obj.app.button(1:3).Visible] = deal(true);
                    [obj.app.button(4:5).Visible] = deal(false);
                    obj.app.tabbutton(1).BackgroundColor = [.94 .94 .94];
                    obj.app.tabbutton(1).FontColor       = 'k';
                    obj.app.tabbutton(2).BackgroundColor = [0 .25 .45];
                    obj.app.tabbutton(2).FontColor       = 'w';
                case 2
                    [obj.app.button(1:3).Visible] = deal(false);
                    [obj.app.button(4:5).Visible] = deal(true);
                    obj.app.tabbutton(1).BackgroundColor = [0 .25 .45];
                    obj.app.tabbutton(1).FontColor       = 'w';
                    obj.app.tabbutton(2).BackgroundColor = [.94 .94 .94];
                    obj.app.tabbutton(2).FontColor       = 'k';
            end
            obj.app.maingrid.RowHeight = {25,80,'1x'};
        end
        %% Fold or unfold the menu
        function cb_fold(obj,~,~)
            if strcmp(obj.app.tip.ImageSource,fullfile(obj.iconFolder,'tip1.png'))
                obj.app.tip.ImageSource = fullfile(obj.iconFolder,'tip2.png');
            else
                obj.app.tip.ImageSource = fullfile(obj.iconFolder,'tip1.png');
                obj.app.tabbutton(1).BackgroundColor = [0 .25 .45];
                obj.app.tabbutton(1).FontColor       = 'w';
                obj.app.tabbutton(2).BackgroundColor = [0 .25 .45];
                obj.app.tabbutton(2).FontColor       = 'w';
                obj.app.maingrid.RowHeight           = {25,0,'1x'};
            end
        end
        %% Hide the menu when moving out of it
        function cb_motion(obj,~,~)
            if obj.app.maingrid.RowHeight{2} > 0 && strcmp(obj.app.tip.ImageSource,fullfile(obj.iconFolder,'tip1.png')) && obj.app.figure.CurrentPoint(2) < obj.app.figure.Position(4)-105
                obj.app.tabbutton(1).BackgroundColor = [0 .25 .45];
                obj.app.tabbutton(1).FontColor       = 'w';
                obj.app.tabbutton(2).BackgroundColor = [0 .25 .45];
                obj.app.tabbutton(2).FontColor       = 'w';
                obj.app.maingrid.RowHeight           = {25,0,'1x'};
            end
        end
        %% Read the function lists
        function readList(obj)
            % Read the algorithm list
            LabelStr    = {'none','single','multi','many','real','binary','permutation','large','constrained','expensive','multimodal','sparse','preference'};
            obj.algList = obj.readList2('Algorithms',LabelStr);
            obj.proList = obj.readList2('Problems',LabelStr);
            obj.metList = obj.readList2('Metrics',{'min','max'});
        end
        function List = readList2(obj,folder,LabelStr)
            List    = {};
            Folders = split(genpath(fullfile(fileparts(mfilename('fullpath')),'..',folder)),pathsep);
            for i = 1 : length(Folders)
                Files = what(Folders{i});
                Files = Files.m;
                for j = 1 : length(Files)
                    f = fopen(Files{j});
                    fgetl(f);
                    str = regexprep(fgetl(f),'^\s*%\s*','','once');
                    fclose(f);
                    labelstr = regexp(str,'(?<=<).*?(?=>)','match');
                    if ~isempty(labelstr)
                        label = false(length(labelstr),length(LabelStr));
                        for k = 1 : length(labelstr)
                            label(k,:) = ismember(LabelStr,split(labelstr{k},'/'));
                        end
                        if any(label(:))
                            List = [List;{label},Files{j}(1:end-2)];
                        end
                    end
                end
            end
        end
        %% Change the module
        function cb_module(obj,~,~,index)
            name = {'module_test','module_app','module_exp'};
            for i = 1 : length(name)
                if isfield(obj.app,name{i})
                    obj.app.(name{i}).app.maingrid.Visible = i==index;
                elseif i == index
                    obj.app.(name{i}) = feval(name{i},obj);
                end
            end
        end
        %% Show the figure of about PlatEMO
        function cb_author(obj,~,~)
            P = obj.app.figure.Position;
            C = obj.app.figure.CurrentPoint;
            f = uifigure('Name','About PlatEMO','Position',[P(1)+C(1) P(2)+C(2)-220 290 200],'Color','w','icon',fullfile(obj.iconFolder,'logo1.png'),'Resize','off');
            uiimage(f,'Position',[145 40 150 150],'ImageSource',fullfile(obj.iconFolder,'logo2.png'));
            uilabel(f,'Position',[10 160 170 30],'Text',obj.app.figure.Name,'HorizontalAlignment','left','FontSize',18);
            uilabel(f,'Position',[10 130 160 40],'Text',{'Evolutionary Multi-Objective','Optimization Platform'},'HorizontalAlignment','left','FontSize',12,'FontAngle','italic','FontColor',[.4 .4 .4]);
            uilabel(f,'Position',[10 80 150 40],'Text',{'Copyright (c) 2021','BIMK Group'},'HorizontalAlignment','left','FontSize',12);
            uilabel(f,'Position',[10 20 150 50],'Text',{'<Contact>','field910921@gmail.com','bimk.ahu.edu.cn'},'HorizontalAlignment','left','FontSize',12);
            uibutton(f,'Position',[210 10 70 22],'Text','OK','ButtonpushedFcn',@(~,~)delete(f));
        end
    end
    methods(Static)
        %% Generate a component
        function app = APP(row,column,app)
            app.Layout.Row    = row;
            app.Layout.Column = column;
        end
    end
end