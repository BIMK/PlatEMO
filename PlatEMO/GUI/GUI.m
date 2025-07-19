classdef GUI < handle
%GUI - The class of the main figure of PlatEMO.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        app = struct();     % All the components
        icon;               % Icons
        algList;            % Algorithm list
        proList;            % Problem list
        metList;            % Metric list
    end
    methods
        %% Establish the figure window
        function obj = GUI()
            % Load the data
            obj.icon = load(fullfile(fileparts(mfilename('fullpath')),'GUI'),'-mat');
            obj.readList();
            
            % Create the window
            obj.app.figure   = uifigure('Name','PlatEMO v4.13','Position',[0 0 1200 650],'Interruptible','off','icon',obj.icon.logo1,'BusyAction','cancel','Visible','off','WindowButtonMotionFcn',@(~,~)[]);
            obj.app.maingrid = uigridlayout(obj.app.figure,'RowHeight',{25,80,'1x'},'ColumnWidth',{'1x'},'Padding',[0 0 0 0],'RowSpacing',0);
            
            % Create the tab buttons
            obj.app.grid(1)    = GUI.APP(1,1,uigridlayout(obj.app.maingrid,'RowHeight',{'1x'},'ColumnWidth',{80,80,'1x',100},'Padding',[0 0 0 0],'ColumnSpacing',0,'BackgroundColor',[0 .25 .45]));
            tempPanel          = GUI.APP(1,1,uipanel(obj.app.grid(1),'BorderType','none','BackgroundColor',[0 .25 .45]));
            obj.app.buttonT(1) = uibutton(tempPanel,'Position',[-5 -5 90 35],'Text','Modules','FontSize',14,'FontColor','k','BackgroundColor',[.94 .94 .94],'ButtonpushedFcn',{@obj.cb_tab,1});
            tempPanel          = GUI.APP(1,2,uipanel(obj.app.grid(1),'BorderType','none','BackgroundColor',[0 .25 .45]));
            obj.app.buttonT(2) = uibutton(tempPanel,'Position',[-5 -5 90 35],'Text','Support','FontSize',14,'FontColor','w','BackgroundColor',[0 .25 .45],'ButtonpushedFcn',{@obj.cb_tab,2});
            tempPanel          = GUI.APP(1,4,uipanel(obj.app.grid(1),'BorderType','none','BackgroundColor',[0 .25 .45],'AutoResizeChildren',false));
            tempImage          = uiimage(tempPanel,'Position',[0 3 99 21],'ImageSource',obj.icon.bar,'ScaleMethod','fill');
            tempPanel2         = uipanel(tempPanel,'Position',[13 4 18 18],'BorderType','none','BackgroundColor',[.549 .6627 .7529]);
            obj.app.buttonT(3) = uibutton(tempPanel2,'Position',[-2 -2 24 24],'Text','','Icon',obj.icon.researchgate,'BackgroundColor',[.549 .6627 .7529],'Tooltip','ResearchGate','ButtonpushedFcn',@(~,~)web('https://www.researchgate.net/publication/365027806','-browser'));
            tempPanel2         = uipanel(tempPanel,'Position',[43 4 18 18],'BorderType','none','BackgroundColor',[.549 .6627 .7529]);
            obj.app.buttonT(4) = uibutton(tempPanel2,'Position',[-2 -2 24 24],'Text','','Icon',obj.icon.github,'BackgroundColor',[.549 .6627 .7529],'Tooltip','GitHub','ButtonpushedFcn',@(~,~)web('https://github.com/BIMK/PlatEMO','-browser'));
            tempPanel2         = uipanel(tempPanel,'Position',[73 4 18 18],'BorderType','none','BackgroundColor',[.549 .6627 .7529]);
            obj.app.buttonT(5) = uibutton(tempPanel2,'Position',[-2 -2 24 24],'Text','','Icon',obj.icon.qq,'BackgroundColor',[.549 .6627 .7529],'Tooltip','QQ','ButtonpushedFcn',@(~,~)web('https://qm.qq.com/q/CmVRo5Adyw','-browser'));
            
            % Create the menu
            obj.app.grid(2)   = GUI.APP(2,1,uigridlayout(obj.app.maingrid,'RowHeight',{'1x',13,1},'ColumnWidth',{1,75,75,75,75,'1x',250,13,1},'Padding',[0 0 0 5],'RowSpacing',5));
            obj.app.button(1) = GUI.APP([1 2],2,uibutton(obj.app.grid(2),'Text',{'Test','Module'},'VerticalAlignment','bottom','FontSize',11,'Icon',obj.icon.test,'IconAlignment','top','Tooltip',{'Test one algorithm on a problem with specified parameter settings.','You can analyse the result and study the performance of the algorithm from various aspects.'},'ButtonpushedFcn',{@obj.cb_module,1}));
            obj.app.button(2) = GUI.APP([1 2],3,uibutton(obj.app.grid(2),'Text',{'Application','Module'},'VerticalAlignment','bottom','FontSize',11,'Icon',obj.icon.application,'IconAlignment','top','Tooltip',{'Use algorithms to solve your own problem.','You can design your own problem and solve it by the suggested algorithms.'},'ButtonpushedFcn',{@obj.cb_module,2}));
            obj.app.button(3) = GUI.APP([1 2],4,uibutton(obj.app.grid(2),'Text',{'Experiment','Module'},'VerticalAlignment','bottom','FontSize',11,'Icon',obj.icon.experiment,'IconAlignment','top','Tooltip',{'Do experiment on multiple algorithms and problems.','You can observe the statistical results shown in a table and save it as an Excel or LaTeX table.'},'ButtonpushedFcn',{@obj.cb_module,3}));
            obj.app.button(4) = GUI.APP([1 2],5,uibutton(obj.app.grid(2),'Text',{'Creation','Module'},'VerticalAlignment','bottom','FontSize',11,'Icon',obj.icon.creation,'IconAlignment','top','Tooltip',{'Create new algorithms via blocks.','You can visually create a new algorithm by connecting blocks and train it on problems.'},'ButtonpushedFcn',{@obj.cb_module,4}));
            obj.app.button(5) = GUI.APP([1 2],2,uibutton(obj.app.grid(2),'Visible',false,'Text',{'About','PlatEMO'},'VerticalAlignment','bottom','FontSize',11,'Icon',obj.icon.author,'IconAlignment','top','ButtonpushedFcn',@obj.cb_author));
            obj.app.button(6) = GUI.APP([1 2],3,uibutton(obj.app.grid(2),'Visible',false,'Text',{'User','Manual'},'VerticalAlignment','bottom','FontSize',11,'Icon',obj.icon.help,'IconAlignment','top','ButtonpushedFcn',@(~,~)web(['file://',fullfile(fileparts(fileparts(mfilename('fullpath'))),'manual.pdf')],'-browser')));
            obj.app.tip       = GUI.APP(2,8,uiimage(obj.app.grid(2),'ImageSource',obj.icon.tip2,'ImageClickedFcn',@obj.cb_fold,'UserData',true));
            tempLine          = GUI.APP(3,[1 9],uipanel(obj.app.grid(2),'BackgroundColor',[.8 .8 .8]));
            
            % Create the modules
            movegui(obj.app.figure,'center');
            obj.app.figure.addlistener('CurrentPoint','PostSet',@obj.cb_motion);
            obj.cb_module([],[],obj.icon.GUIsetting);
            obj.app.figure.Visible = 'on';
            
            % Show images
            index = num2str(randi(3));
            if isfield(obj.icon,['image',index])
                GUI.APP([1 2],7,uiimage(obj.app.grid(2),'ImageSource',obj.icon.(['image',index]),'ImageClickedFcn',@(~,~)web(['https://bimk.github.io/Conference-Competition/?page=',index],'-browser')));
            end
        end
    end
	methods(Access = private)
        %% Change the menu buttons
        function cb_tab(obj,~,~,type)
            switch type
                case 1
                    [obj.app.button(1:4).Visible] = deal(true);
                    [obj.app.button(5:6).Visible] = deal(false);
                    obj.app.buttonT(1).BackgroundColor = [.94 .94 .94];
                    obj.app.buttonT(1).FontColor       = 'k';
                    obj.app.buttonT(2).BackgroundColor = [0 .25 .45];
                    obj.app.buttonT(2).FontColor       = 'w';
                case 2
                    [obj.app.button(1:4).Visible] = deal(false);
                    [obj.app.button(5:6).Visible] = deal(true);
                    obj.app.buttonT(1).BackgroundColor = [0 .25 .45];
                    obj.app.buttonT(1).FontColor       = 'w';
                    obj.app.buttonT(2).BackgroundColor = [.94 .94 .94];
                    obj.app.buttonT(2).FontColor       = 'k';
            end
            obj.app.maingrid.RowHeight = {25,80,'1x'};
        end
        %% Fold or unfold the menu
        function cb_fold(obj,~,~)
            obj.app.tip.UserData = ~obj.app.tip.UserData;
            if obj.app.tip.UserData
                obj.app.tip.ImageSource = obj.icon.tip2;
            else
                obj.app.tip.ImageSource = obj.icon.tip1;
                obj.app.buttonT(1).BackgroundColor = [0 .25 .45];
                obj.app.buttonT(1).FontColor       = 'w';
                obj.app.buttonT(2).BackgroundColor = [0 .25 .45];
                obj.app.buttonT(2).FontColor       = 'w';
                obj.app.maingrid.RowHeight         = {25,0,'1x'};
            end
        end
        %% Hide the menu when moving out of it
        function cb_motion(obj,~,~)
            if obj.app.maingrid.RowHeight{2} > 0 && ~obj.app.tip.UserData && obj.app.figure.CurrentPoint(2) < obj.app.figure.Position(4)-105
                obj.app.buttonT(1).BackgroundColor = [0 .25 .45];
                obj.app.buttonT(1).FontColor       = 'w';
                obj.app.buttonT(2).BackgroundColor = [0 .25 .45];
                obj.app.buttonT(2).FontColor       = 'w';
                obj.app.maingrid.RowHeight         = {25,0,'1x'};
            end
        end
        %% Read the function lists
        function readList(obj)
            LabelStr    = {'none','single','multi','many','real','integer','label','binary','permutation','large','constrained','expensive','multimodal','sparse','dynamic','multitask','bilevel','robust'};
            obj.algList = obj.readList2('Algorithms',LabelStr);
            obj.proList = obj.readList2('Problems',LabelStr);
            obj.metList = obj.readList2('Metrics',[LabelStr,'min','max']);
        end
        function List = readList2(obj,folder,LabelStr)
            List    = {};
            Folders = split(genpath(fullfile(fileparts(mfilename('fullpath')),'..',folder)),pathsep);
            for i = 1 : length(Folders) - 1 
                Files = what(Folders{i});
                Files = Files.m;
                for j = 1 : length(Files)
                    try
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
                                emptyLine = find(~any(label,2));
                                if isempty(emptyLine)
                                    year = ' ';
                                else
                                    year = labelstr{emptyLine(1)};
                                    label(emptyLine,:) = [];
                                end
                                List = [List;{label},Files{j}(1:end-2),year];
                            end
                        end
                    catch
                    end
                end
            end
        end
        %% Change the module
        function cb_module(obj,~,~,GUIsetting)
            name = {'module_test','module_app','module_exp','module_cre'};
            for i = 1 : length(name)
                if isfield(obj.app,name{i})
                    obj.app.(name{i}).app.maingrid.Visible = i==GUIsetting;
                elseif i == GUIsetting
                    obj.app.(name{i}) = feval(name{i},obj);
                end
            end
            save(fullfile(fileparts(mfilename('fullpath')),'GUI'),'GUIsetting','-append');
        end
        %% Show the figure of about PlatEMO
        function cb_author(obj,~,~)
            P = obj.app.figure.Position;
            C = obj.app.figure.CurrentPoint;
            f = uifigure('Name','About PlatEMO','Position',[P(1)+C(1) P(2)+C(2)-220 290 200],'Color','w','icon',obj.icon.logo1,'WindowStyle','modal','Resize','off');
            uiimage(f,'Position',[145 40 150 150],'ImageSource',obj.icon.logo2);
            uilabel(f,'Position',[10 160 170 30],'Text',obj.app.figure.Name,'HorizontalAlignment','left','FontSize',18);
            uilabel(f,'Position',[10 130 160 40],'Text',{'Evolutionary Multi-Objective','Optimization Platform'},'HorizontalAlignment','left','FontSize',12,'FontAngle','italic','FontColor',[.4 .4 .4]);
            uilabel(f,'Position',[10 80 150 40],'Text',{'Copyright (c) 2024','BIMK Group'},'HorizontalAlignment','left','FontSize',12);
            uilabel(f,'Position',[10 20 150 50],'Text',{'<Contact>','field910921@gmail.com','bimk.ahu.edu.cn'},'HorizontalAlignment','left','FontSize',12);
            uibutton(f,'Position',[210 10 70 22],'Text','OK','ButtonpushedFcn',@(~,~)delete(f));
        end
        %% Load images after closing the GUI
        function delete(obj)
            imagetime = clock;
            if ~isfield(obj.icon,'imagetime') || etime(imagetime,obj.icon.imagetime) > 259200
                try
                    image1 = webread('https://bimk.github.io/imgs/picture1.jpg');
                    image2 = webread('https://bimk.github.io/imgs/picture2.jpg');
                    image3 = webread('https://bimk.github.io/imgs/picture3.jpg');
                    save(fullfile(fileparts(mfilename('fullpath')),'GUI'),'image1','-append');
                    save(fullfile(fileparts(mfilename('fullpath')),'GUI'),'image2','-append');
                    save(fullfile(fileparts(mfilename('fullpath')),'GUI'),'image3','-append');
                    save(fullfile(fileparts(mfilename('fullpath')),'GUI'),'imagetime','-append');
                catch
                end
            end
        end
    end
    methods(Static)
        %% Generate a component
        function app = APP(row,column,app)
            app.Layout.Row    = row;
            app.Layout.Column = column;
        end
        %% Generate state buttons of labels
        function [stateButton,label] = GenerateLabelButton(grid,values,cbFcn)
            label(1) = GUI.APP(1,[1 3],uilabel(grid,'Text','Number of objectives','VerticalAlignment','bottom','FontSize',12,'FontColor',[.15 .6 .2],'FontWeight','bold'));
            stateButton(1)  = GUI.APP(2,1,uibutton(grid,'state','Text','single','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',values(1),'Tooltip','The problem has a single objective','ValueChangedFcn',{cbFcn,1}));
            stateButton(2)  = GUI.APP(2,2,uibutton(grid,'state','Text','multi','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',values(2),'Tooltip','The problem has 2 or 3 objectives','ValueChangedFcn',{cbFcn,2}));
            stateButton(3)  = GUI.APP(2,3,uibutton(grid,'state','Text','many','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',values(3),'Tooltip','The problem has more than 3 objectives','ValueChangedFcn',{cbFcn,3}));
            label(2) = GUI.APP(3,[1 3],uilabel(grid,'Text','Encoding scheme','VerticalAlignment','bottom','FontSize',12,'FontColor',[.15 .6 .2],'FontWeight','bold'));
            stateButton(4)  = GUI.APP(4,1,uibutton(grid,'state','Text','real','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',values(4),'Tooltip','The decision variables are real numbers','ValueChangedFcn',{cbFcn,4}));
            stateButton(5)  = GUI.APP(4,2,uibutton(grid,'state','Text','integer','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',values(5),'Tooltip','The decision variables are integers','ValueChangedFcn',{cbFcn,5}));
            stateButton(6)  = GUI.APP(4,3,uibutton(grid,'state','Text','label','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',values(6),'Tooltip','The decision variables are labels','ValueChangedFcn',{cbFcn,6}));
            stateButton(7)  = GUI.APP(5,1,uibutton(grid,'state','Text','binary','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',values(7),'Tooltip','The decision variables are binary numbers','ValueChangedFcn',{cbFcn,7}));
            stateButton(8)  = GUI.APP(5,2,uibutton(grid,'state','Text','permutation','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',values(8),'Tooltip','The decision vector is a permutation','ValueChangedFcn',{cbFcn,8}));
            label(3) = GUI.APP(6,[1 3],uilabel(grid,'Text','Special difficulties','VerticalAlignment','bottom','FontSize',12,'FontColor',[.15 .6 .2],'FontWeight','bold'));
            stateButton(9)  = GUI.APP(7,1,uibutton(grid,'state','Text','large','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',values(9),'Tooltip','The problem has more than 100 decision variables','ValueChangedFcn',{cbFcn,9}));
            stateButton(10) = GUI.APP(7,2,uibutton(grid,'state','Text','constrained','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',values(10),'Tooltip','The problem has constraints','ValueChangedFcn',{cbFcn,10}));
            stateButton(11) = GUI.APP(7,3,uibutton(grid,'state','Text','expensive','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',values(11),'Tooltip','The objectives are computationally time-consuming','ValueChangedFcn',{cbFcn,11}));
            stateButton(12) = GUI.APP(8,1,uibutton(grid,'state','Text','multimodal','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',values(12),'Tooltip','The objectives are multimodal','ValueChangedFcn',{cbFcn,12}));
            stateButton(13) = GUI.APP(8,2,uibutton(grid,'state','Text','sparse','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',values(13),'Tooltip','Most decision variables of the optimal solutions are zero','ValueChangedFcn',{cbFcn,13}));
            stateButton(14) = GUI.APP(8,3,uibutton(grid,'state','Text','dynamic','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',values(14),'Tooltip','The objectives vary periodically','ValueChangedFcn',{cbFcn,14}));
            stateButton(15) = GUI.APP(9,1,uibutton(grid,'state','Text','multitask','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',values(15),'Tooltip','The problem has multiple tasks to be solved simultaneously','ValueChangedFcn',{cbFcn,15}));
            stateButton(16) = GUI.APP(9,2,uibutton(grid,'state','Text','bilevel','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',values(16),'Tooltip','The problem has two nested objectives','ValueChangedFcn',{cbFcn,16}));
            stateButton(17) = GUI.APP(9,3,uibutton(grid,'state','Text','robust','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',values(17),'Tooltip','The objectives are influenced by uncertain factors','ValueChangedFcn',{cbFcn,17}));
        end
        %% Update the list of algorithms and problems
        function func = UpdateAlgProList(index,stateButton,varargin)
            if index > 0
                if index < 4
                    [stateButton(1:3).Value] = deal(0);
                    stateButton(index).Value = 1;
                elseif index < 9
                    [stateButton(4:8).Value] = deal(0);
                    stateButton(index).Value = 1;
                end
            end
            filter = [stateButton.Value];
            func   = @(s)all(any(repmat([true,filter],size(s,1),1)&s,2)) && all((any(s(:,2:end),1)&filter)==filter);
            for i = 1 : 4 : length(varargin)
                [~,drop,~,list] = deal(varargin{i:i+3});
                drop.UserData   = find(cellfun(func,list(:,1)));
                drop.Items      = flip([unique(list(drop.UserData,3));'All year']);
            end
            GUI.UpdateAlgProListYear(varargin{:});
        end
        %% Update the list of algorithms and problems by year
        function UpdateAlgProListYear(varargin)
            for i = 1 : 4 : length(varargin)
                [listBox,drop,label,list] = deal(varargin{i:i+3});
                if strcmp(drop.Value,'All year')
                    listBox.Items = ['(Open File)';list(drop.UserData,2)];
                    listBox.Value = {};
                    label.Text    = sprintf('%d / %d',length(drop.UserData),size(list,1));
                else
                    index = ismember(list(drop.UserData,3),drop.Value);
                    listBox.Items = ['(Open File)';list(drop.UserData(index),2)];
                    listBox.Value = {};
                    label.Text    = sprintf('%d / %d',length(drop.UserData(index)),size(list,1));
                end
            end
        end
        %% Update the parameter list of algorithms and problems
        function UpdateAlgProPara(fig,list,paraList,fileType,paraType)
            filename = list.Value;
            if contains(filename,'Open File')
                [Name,Path] = uigetfile({'*.m','MATLAB class'});
                figure(fig);
                if Name ~= 0
                    try
                        filename = fullfile(Path,Name);
                        f   = fopen(filename);
                        str = fgetl(f);
                        fclose(f);
                        assert(contains(str,['< ',fileType]));
                        addpath(Path);
                    catch
                        uialert(fig,['The selected file is not a subclass of ',fileType,'.'],'Error');
                        return;
                    end
                else
                    return;
                end
            else
                filename = [filename,'.m'];
            end
            if abs(paraType) ~= 2
                paraList.del([],paraType);
            end
            paraList.add(filename,paraType);
            paraList.flush();
        end
        %% Read parameter settings from the parameter list
        function [name,para] = GetParameterSetting(listItem)
            name = listItem.title.Text;
            para = cell(1,length(listItem.edit));
            for j = 1 : length(para)
                if ~isempty(listItem.edit(j).Value)
                    if listItem.edit(j).Value(1) == '''' || listItem.edit(j).Value(1) == '"'
                        para{j} = listItem.edit(j).Value(2:end-1);
                    else
                        para{j} = str2num(listItem.edit(j).Value);
                    end
                    assert(~isempty(para{j}),'the parameter "%s" of %s is illegal.',listItem.label(j).Text,listItem.title.Text);
                end
            end
        end
        %% Save a population
        function SavePopulation(fig,Population,type)
            if type == 1
                Population = Population.best;
            end
            if isempty(Population)
                uialert(fig,'No solution can be saved, since all solutions are infeasible.','Error');
            else
                Data = [Population.decs,Population.objs,Population.cons];
                try
                    [Name,Path] = uiputfile({'*.txt','Text file';'*.dat','Text file';'*.csv','Text file';'*.mat','MAT file';'*.xlsx','Excel table'},'','data');
                    figure(fig);
                    if ischar(Name)
                        [~,~,Type] = fileparts(Name);
                        switch Type
                            case '.mat'
                                save(fullfile(Path,Name),'Data','-mat');
                            otherwise
                                writematrix(Data,fullfile(Path,Name));
                        end
                    end
                catch err
                    uialert(fig,'Fail to save the result, please refer to the command window for details.','Error');
                    rethrow(err);
                end
            end
        end
    end
end