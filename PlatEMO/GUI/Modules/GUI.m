classdef GUI < handle
%GUI - The class of the main figure window of the platform.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private, SetObservable)
        figure;             % The main figure
        menu;               % The main menu
        control = struct(); % All the controls
        icons;              % All the icons
        algList;            % Algorithm list
        proList;            % Problem list
        metList;            % Metric list
    end
    methods
        %% Establish the figure window
        function obj = GUI()
            % Read the icons
            load(fullfile('GUI','Modules','data'),'-mat');
            obj.icons = icons;
            % Read the function lists
            obj.readList();
            % Create the figure
            obj.figure = newFigure([0 0 1200 650],'PlatEMO v2.0');
            % Create the menu
            obj.menu = newTab(obj.figure,[0 551 1202 100,1 1 0 1],{'Modules','Help'});
            obj.figure.busy = false;
            % The first tab of the menu
            obj.control.summaryLabel    = newLabel(obj.menu.panels(1),[700 3 500 65,0 1 1 0],'','FontSize',8,'HorizontalAlignment','left');
            obj.control.moduleButton(1) = newButtonC(obj.menu.panels(1),[15 3 75 65,1 0 1 0],obj.icons.test,{'','','Test','module'},'choosed',true,'UserData',@module_test,...
                                                     'callback',@obj.changeModule,'movecallback',@obj.changeSummary,'moveoutcallback',@obj.changeSummary);
            obj.control.moduleButton(2) = newButtonC(obj.menu.panels(1),[95 3 75 65,1 0 1 0],obj.icons.experiment,{'','','Experiment','module'},'choosed',true,'UserData',@module_experiment,...
                                                     'callback',@obj.changeModule,'movecallback',@obj.changeSummary,'moveoutcallback',@obj.changeSummary);
            obj.changeModule(obj.control.moduleButton(1));
            % The second tab of the menu
            newButtonC(obj.menu.panels(2),[15 3 50 65,1 0 1 0],obj.icons.author,{'','','About','PlatEMO'},'callback',@obj.cb_author);
            newLine(obj.menu.panels(2),[75 5 1 60,1 0 1 1]);
            newButtonC(obj.menu.panels(2),[85 3 50 65,1 0 1 0],obj.icons.help,{'','','How to','use'},'callback',@obj.cb_help);
        end
    end
	methods(Access = private)
        %% Read the function lists
        function readList(obj)
            obj.algList = obj.readList2('Algorithms','algorithm');
            obj.proList = obj.readList2('Problems','problem');
            obj.metList = obj.readList2('Metrics','metric');
            obj.algList = [{'Highlighted'},{{'ARMOEA','GrEA','IBEA','MOEAD','NSGAII','NSGAIISDR','NSGAIII','RVEA'}};obj.algList];
        end
        function List = readList2(obj,folder,class)
            Folders  = regexp(genpath(fullfile(folder)),'[;,:]','split');
            FuncList = {};
            for i = 1 : length(Folders)-1
                Files = what(Folders{i});
                Files = Files.m;
                for j = 1 : length(Files)
                    f = fopen(Files{j});
                    fgetl(f);
                    str = regexprep(fgetl(f),'^\s*%\s*','','once');
                    fclose(f);
                    labels = regexp(str,'(?<=<).*?(?=>)','match');
                    if ~isempty(labels) && strcmp(labels{1},class)
                        if length(labels) == 1
                            labels{2} = '';
                        end
                        FuncList = [FuncList;labels(2),Files(j)];
                    end
                end
            end
            FuncList(:,2) = cellfun(@(S)S(1:end-2),FuncList(:,2),'UniformOutput',false);
            ULabel        = unique(FuncList(:,1));
            [~,Index]     = ismember(FuncList(:,1),ULabel);
            List          = cell(length(ULabel),2);
            for i = 1 : length(ULabel)
                List(i,:) = [ULabel(i),{FuncList(Index==i,2)}];
            end
        end
        %% Change the current module
        function changeModule(obj,hObject,eventdata)
            for h = obj.control.moduleButton
                if ~isequal(h,hObject)
                    h.value = false;
                    if isa(h.handle.UserData,'module')
                        h.handle.UserData.panel.visible = false;
                    end
                else
                    h.value = true;
                    if isa(h.handle.UserData,'module')
                        h.handle.UserData.panel.visible = true;
                    else
                        h.handle.UserData = h.handle.UserData(obj);
                        h.handle.UserData.panel.position(3:4) = h.handle.UserData.panel.position(3:4) + obj.figure.handle.Position(3:4) - [1200 650];
                    end
                end
            end
        end
        %% Change the summary
        function changeSummary(obj,hObject,eventdata)
            current = find([obj.control.moduleButton.moved],1);
            if ~isempty(current)
                Data = obj.control.moduleButton(current).handle.UserData;
                if isa(Data,'module')
                    obj.control.summaryLabel.handle.String = Data.summary;
                else
                    obj.control.summaryLabel.handle.String = eval([func2str(Data),'.summary']);
                end
            else
                obj.control.summaryLabel.handle.String = '';
            end
        end
        %% Show the authors
        function cb_author(obj,hObject,eventdata)
            sf = newFigure([1 1 300 220],'About PlatEMO',obj.figure);
            sf.handle.Color = [.95 .95 1];
            newIcon(sf,[155 60 150 150,1 0 1 0],obj.icons.logo);
            newLabel(sf,[10 180 170 30,1 0 1 0],obj.figure.handle.Name,'HorizontalAlignment','left','FontSize',13);
            newLabel(sf,[10 150 160 40,1 0 1 0],'Evolutionary Multi-Objective Optimization Platform','HorizontalAlignment','left','FontSize',8,'FontAngle','italic','ForegroundColor',[.4 .4 .4]);
            newLabel(sf,[10 110 150 40,1 0 1 0],'Copyright (c) 2018-2019 BIMK Group','HorizontalAlignment','left','FontSize',9);
            newLabel(sf,[10 90 150 15,1 0 1 0],'<Contact>','HorizontalAlignment','left','FontSize',8);
            newLink(sf,[10 70 120 20,1 0 1 0],'field910921@gmail.com','FontSize',11);
            newLink(sf,[10 50 120 20,1 0 1 0],'bimk.ahu.edu.cn','FontSize',11);
            newLabel(sf,[10 30 150 15,1 0 1 0],'<Join us on GitHub>','HorizontalAlignment','left','FontSize',8);
            newLink(sf,[10 10 150 20,1 0 1 0],'github.com/BIMK/PlatEMO','FontSize',11);
            newButton(sf,[220 15 70 22,1 0 1 0],'OK','callback',@sf.WindowClose);
            sf.busy = false;
        end
        %% Show the guidance
        function cb_help(obj,hObject,eventdata)
            current = find([obj.control.moduleButton.value],1);
            obj.control.moduleButton(current).handle.UserData.showGuidance();
        end
    end
end