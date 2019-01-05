classdef module_experiment < module
%module_experiment - Experiment module.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods(Access = ?GUI)
        %% Constructor
        function obj = module_experiment(GUI)
            Panel = newMpanel(GUI.figure,[1 1 1200 549,1 1 1 1],3);
            Panel.move(1,-246);
            Panel.move(2,-442);
            obj = obj@module(GUI,Panel);
            
            % The first panel
            panel = Panel.panels(1);
            note = newLabel(panel,[10 10 140 115,0 0 1 1],'','FontSize',9,'HorizontalAlignment','left');
            newLabel(panel,[15 530 135 15,0 0 0 1],'Function selection','ForegroundColor',[.7 .7 .7],'FontSize',9);
            obj.control.setLabel(1)  = newLabelButton(panel,[10 460 135 60,0 0 0 1],'Algorithm','ForegroundColor',[.2 .4 .7],'FontSize',12,'FontWeight','bold','choosed',true,...
                                                      'callback',@(~,~)obj.cb_setLabel(1),'movecallback',@(~,~)set(note.handle,'String','Select the algorithms to be executed'));
            obj.control.setLabel(2)  = newLabelButton(panel,[10 395 135 60,0 0 0 1],'Problem','ForegroundColor',[.9 .5 .2],'FontSize',12,'FontWeight','bold','choosed',true,...
                                                      'callback',@(~,~)obj.cb_setLabel(2),'movecallback',@(~,~)set(note.handle,'String','Select the problems to be tested'));
            obj.control.algPopmenu   = newPopmenu(panel,[20 470 115 20,0 0 0 1],[.2 .4 .7],obj.GUI.algList,'callback',@(~,~)obj.cb_setList('alg'));
            obj.control.proPopmenu   = newPopmenu(panel,[20 405 115 20,0 0 0 1],[.9 .5 .2],obj.GUI.proList,'callback',@(~,~)obj.cb_setList('pro'));
            obj.control.numpopLabel  = newLabel(panel,[10 335 135 47,0 0 0 1],'Number of results','ForegroundColor',[.3 .5 .2],'FontSize',10,'FontWeight','bold',...
                                                      'movecallback',@(~,~)set(note.handle,'String','Set the number of populations saved during the evolution of each run'));
            obj.control.runtimeLabel = newLabel(panel,[10 275 135 47,0 0 0 1],'Number of runs','ForegroundColor',[.3 .5 .2],'FontSize',11,'FontWeight','bold',...
                                                      'movecallback',@(~,~)set(note.handle,'String','Set the number of runs'));
            obj.control.folderLabel  = newLabel(panel,[10 210 135 47,0 0 0 1],'File path','ForegroundColor',[.3 .5 .2],'FontSize',11,'FontWeight','bold',...
                                                      'movecallback',@(~,~)set(note.handle,'String','Set the file for saving experimental settings, all the results are saved in the same folder'));
            obj.control.numpopEdit   = newEdit(panel,[20 335 115 20,0 0 0 1],'String','1');
            obj.control.runtimeEdit  = newEdit(panel,[20 275 115 20,0 0 0 1],'String','30');
            obj.control.folderEdit   = newEdit(panel,[20 210 105 20,0 0 0 1],'String',fullfile('Data','Setting.mat'),'FontSize',9);
            obj.control.folderButton = newButtonT(panel,[128 213 17 17,0 0 0 1],obj.GUI.icons.load,'Load experimental settings','callback',@(~,~)obj.cb_slSetting());
            obj.control.runButton    = newButton(panel,[20 150 50 40,0 0 0 1],'RUN','FontSize',15,'enable',false,'callback',@(~,~)obj.cb_run(0));
            obj.control.runButton2   = newButton(panel,[75 150 60 40,0 0 0 1],'<html><center>RUN in<br>parallel</center></html>','FontSize',11,'enable',false,'callback',@(~,~)obj.cb_run(1));
            newLine(panel,[5 130 145 1,0 0 0 1]);
            
            % The second panel
            panel = Panel.panels(2);
            newLabel(panel,[15 530 180 15,0 0 0 1],'Parameter setting','ForegroundColor',[.7 .7 .7],'FontSize',9);
            obj.control.setPanel(1) = ParameterList(panel,[6 10 194 514,1 1 1 1],note,obj.GUI.icons,'callback',@obj.cb_setList2);
            obj.control.setPanel(2) = ParameterList(panel,[6 10 194 514,1 1 1 1],note,obj.GUI.icons,'callback',@obj.cb_setList2);
            obj.cb_setLabel(1);

            % The third panel
            panel = Panel.panels(3);
            newLabel(panel,[369 530 100 15,0 0 0 1],'Result display','ForegroundColor',[.7 .7 .7],'FontSize',9);
            obj.control.tableToolBar   = newPanel(panel,[20 499 798 25,1 1 0 1],[.95 .95 1]);
            obj.control.tableButton(1) = newButtonT(obj.control.tableToolBar,[5 1 25 25,1 0 1 0],obj.GUI.icons.save,'Save table','callback',@(~,~)obj.saveTable());
            newLine(obj.control.tableToolBar,[35 1 1 25,1 0 1 0]);
            obj.control.tableButton(2) = newButtonT(obj.control.tableToolBar,[40 1 25 25,1 0 1 0],'N','Show population size','choosed',true,'FontSize',12,'callback',@(~,~)obj.updataTableColumn());
            obj.control.tableButton(3) = newButtonT(obj.control.tableToolBar,[66 1 25 25,1 0 1 0],'M','Show number of objectives','choosed',true,'value',true,'FontSize',12,'callback',@(~,~)obj.updataTableColumn());
            obj.control.tableButton(4) = newButtonT(obj.control.tableToolBar,[92 1 25 25,1 0 1 0],'D','Show number of decision variables','choosed',true,'value',true,'FontSize',12,'callback',@(~,~)obj.updataTableColumn());
            obj.control.tableButton(5) = newButtonT(obj.control.tableToolBar,[118 1 25 25,1 0 1 0],'FEs','Show number of evaluations','choosed',true,'FontSize',10,'callback',@(~,~)obj.updataTableColumn());
            newLine(obj.control.tableToolBar,[148 1 1 25,1 0 1 0]);
            obj.control.tableButton(6) = newButtonT(obj.control.tableToolBar,[153 1 25 25,1 0 1 0],obj.GUI.icons.std,'Show standard deviation','choosed',true,'value',true,'callback',@(~,~)obj.flushTable);
            obj.control.tableButton(7) = newButtonT(obj.control.tableToolBar,[179 1 25 25,1 0 1 0],obj.GUI.icons.ranksum,'Show ranksum test','choosed',true,'value',true,'callback',@(~,~)obj.flushTable);
            obj.control.ranksumMenu    = newMenu(panel,[1 1 145 20]);
            newTip(obj.control.tableToolBar,[1 1 12 25],obj.control.tableButton(7),'callback',@(~,~)obj.control.ranksumMenu.show());
            obj.control.showButton     = newPopmenu3(obj.control.tableToolBar,[681 1 105 25,0 1 1 0],[{obj.GUI.icons.runTimes,[],obj.GUI.icons.runtime,[]},cell(1,length(obj.GUI.metList{1,2})+length(obj.GUI.metList{2,2}))],...
                                                     [{'Run times','','runtime',''},obj.GUI.metList{2,2}',obj.GUI.metList{1,2}'],'callback',@(~,~)obj.flushTable);
            obj.control.table          = newTable(panel,[24 90 788 390,1 1 1 1],'Fontname','Times New Roman','Fontsize',9,'ButtonDownFcn',@obj.cb_table);
            obj.control.tableMenu      = newMenu(panel,[1 1 245 20]);
            obj.control.tableMenu.add(obj.GUI.icons.showPF,'PF with median metric value','callback',@(~,~)obj.cb_tableMenu(1));
            obj.control.tableMenu.add([],'PF with best metric value','callback',@(~,~)obj.cb_tableMenu(2));
            obj.control.tableMenu.add([],'');
            obj.control.tableMenu.add(obj.GUI.icons.showPS,'PS with median metric value','callback',@(~,~)obj.cb_tableMenu(3));
            obj.control.tableMenu.add([],'PS with best metric value','callback',@(~,~)obj.cb_tableMenu(4));
            obj.control.tableMenu.add([],'');
            obj.control.tableMenu.add(obj.GUI.icons.showMetric,'Variation of average metric values','callback',@(~,~)obj.cb_tableMenu(5));
            obj.control.playbar        = newBar(panel,[24 45 788 27,1 1 1 0]);
            obj.control.playlabel      = newLabel(panel,[24 7 340 20,1 0 1 0],'','FontSize',9,'HorizontalAlignment','left','ForegroundColor',[.4 .4 .4]);
            obj.control.playStart      = newButtonI(panel,[404 5 30 30,0 0 1 0],obj.GUI.icons.start,'enable',false,'callback',@obj.cb_start);
            obj.control.playPause      = newButtonI(panel,[404 5 30 30,0 0 1 0],obj.GUI.icons.pause,'visible',false,'callback',@obj.cb_pause);
            obj.control.playStop       = newButtonI(panel,[364 5 30 30,0 0 1 0],obj.GUI.icons.stop,'callback',@obj.cb_stop);
            obj.control.waitlabel      = newLabel(panel,[380 260 100 20,0 0 0 0],'Please wait ... ...','visible',false);
        end
    end
    methods(Static)
        %% Obtain the summary
        function words = summary()
            words = ['Perform multiple multi-objective optimization algorithms on multiple ',...
            'multi-objective test problems at the same time. The statistical results are ',...
            'shown in a table, and the difference between the performance of the algorithms is ',...
            'highlighted. The results can be saved as Excel table or LaTeX table.'];
        end
        %% Callback of the figure in guidance mode
        function guidance(obj)
            obj.GUI.figure.enable        = false;
            obj.control.guideLabel.state = true;
            obj.control.guideIndex       = obj.control.guideIndex + 1;
            switch obj.control.guideIndex
                case 1
                    obj.updateGuideLabel(' Select several algorithms and set their parameters ',obj.control.setLabel(1));
                    obj.control.algPopmenu.button.state = true;
                    obj.cb_setLabel(1);
                    obj.control.setPanel(1).state = true;
                case 2
                    obj.updateGuideLabel(' Select several problems and set their parameters ',obj.control.setLabel(2));
                    obj.control.proPopmenu.button.state = true;
                    obj.cb_setLabel(2);
                    obj.control.setPanel(2).state = true;
                case 3
                    obj.updateGuideLabel(' Set the number of populations saved in each file ',obj.control.numpopLabel);
                    obj.control.numpopEdit.state = true;
                case 4
                    obj.updateGuideLabel(' Set the number of runs ',obj.control.runtimeLabel);
                    obj.control.runtimeEdit.state = true;
                case 5
                    obj.updateGuideLabel(' Select the folder for saving results ',obj.control.folderLabel);
                    obj.control.folderEdit.state   = true;
                    obj.control.folderButton.state = true;
                case 6
                    obj.updateGuideLabel(' Run the experiment ',obj.control.runButton);
                    obj.control.runButton2.state = true;
                case 7
                    obj.updateGuideLabel(' Display the results from different aspects ',obj.control.tableToolBar);
                case 8
                    obj.updateGuideLabel(' Save the results as Excel table or LaTeX table ',obj.control.tableButton(1));
                otherwise
                    obj.GUI.figure.busy   = false;
                    obj.GUI.figure.enable = true;
                    obj.cb_setLabel(1);
                    delete(obj.control.guideLabel);
            end
        end
    end
    methods(Access = private)
        %% Select the setting panel
        function cb_setLabel(obj,thisIndex)
            [obj.control.setLabel([1:thisIndex-1,thisIndex+1:end]).value] = deal(false);
            obj.control.setLabel(thisIndex).value = true;
            [obj.control.setPanel([1:thisIndex-1,thisIndex+1:end]).visible] = deal(false);
            obj.control.setPanel(thisIndex).visible = true;
        end
        %% Save or load experimental settings
        function success = cb_slSetting(obj,filename)
            if nargin < 2
                % Load experimental settings
                folder        = fileparts(obj.control.folderEdit.string);
                [file,folder] = uigetfile('*.mat','',folder);
                filename      = fullfile(folder,file);
                if exist(filename,'file') == 2
                    try
                        load(filename,'Setting','Environment','-mat');
                        obj.control.setPanel(1).update('replaceall',Setting{1,1},[.2 .4 .7],2);
                        obj.control.setPanel(2).update('replaceall',Setting{2,1},[.9 .5 .2],-2);
                        for i = 1 : 2
                            Edits = [obj.control.setPanel(i).items.edits];
                            if ~isempty(Edits)
                                set([Edits.handle],{'String'},Setting{i,2});
                            end
                        end
                        obj.control.numpopEdit.handle.String  = Environment{1};
                        obj.control.runtimeEdit.handle.String = Environment{2};
                        obj.control.folderEdit.handle.String  = filename;
                    catch
                        errordlg(sprintf('Fail to load experimental settings from %s, since the file does not contain any experimental setting information.',filename),get(gcf,'Name'),'modal'); beep;
                        success = false;
                        return;
                    end
                end
            else
                % Save experimental settings
                Setting = cell(length(obj.control.setPanel),2);
                for i = 1 : size(Setting,1)
                    Setting{i,1} = {obj.control.setPanel(i).items.name};
                    Edits        = [obj.control.setPanel(i).items.edits];
                    if ~isempty(Edits)
                        Setting{i,2} = get([Edits.handle],'String');
                    end
                    if ~iscell(Setting{i,2})
                        Setting{i,2} = {Setting{i,2}};
                    end
                end
                Environment = {obj.control.numpopEdit.handle.String,obj.control.runtimeEdit.handle.String};
                if isempty(filename)
                    filename = fullfile('Data','Setting.mat');
                end
                try
                    [folder,file] = fileparts(filename);
                    if isempty(file)
                        file = 'Setting';
                    end
                    filename = fullfile(folder,[file,'.mat']);
                    if exist(folder,'dir') ~= 7
                        [~,~] = mkdir(folder);
                    end
                    save(filename,'Setting','Environment','-mat');
                    obj.control.folderEdit.handle.String = filename;
                catch
                    errordlg('The file path is illegal, cannot save the experimental settings.',get(gcf,'Name'),'modal'); beep;
                    success = false;
                    return;
                end
            end
            success = true;
        end
        %% Update the parameter setting lists of algorithm and problem
        function cb_setList(obj,type)
            switch type
                case 'alg'
                    obj.control.setPanel(1).update('add',obj.control.algPopmenu.string,[.2 .4 .7],2)
                    obj.cb_setLabel(1);
                case 'pro'
                    obj.control.setPanel(2).update('add',obj.control.proPopmenu.string,[.9 .5 .2],-2);
                    obj.cb_setLabel(2);
            end
        end
        %% Update the state of run button
        function cb_setList2(obj,hObject,eventdata)
            obj.control.runButton.enable  = ~isempty(obj.control.setPanel(1).items) && ~isempty(obj.control.setPanel(2).items);
            obj.control.runButton2.enable = obj.control.runButton.enable;
            obj.control.playStart.enable  = obj.control.runButton.enable;
        end
        %% Run the algorithms
        function cb_run(obj,inparallel)
            % Read all the parameters
            items     = [obj.control.setPanel.items];
            labels    = [items.labels];
            edits     = [items.edits];
            ParaName  = get([labels.handle],'String');
            Parameter = strtrim(get([edits.handle],'String'));
            % Calculate the values of all parameters
            for i = 1 : length(Parameter)
                if ~isempty(Parameter{i})
                    Parameter{i} = str2num(Parameter{i});
                    if isempty(Parameter{i})
                        errordlg(sprintf('The value of parameter <%s> is illegal.',ParaName{i}),get(gcf,'Name'),'modal'); beep;
                        return;
                    end
                end
            end
            % Generate the parameter settings of algorithms and problems
            len1 = cumsum(arrayfun(@(S)length(S.items),obj.control.setPanel));
            len2 = cumsum(arrayfun(@(S)length(S.edits),items));
            len1 = [0,len1];
            len2 = [0,len2];
            parameters = cell(length(len2)-1,1);
            for i = 1 : length(parameters)
                parameters{i} = Parameter(len2(i)+1:len2(i+1))';
            end
            names = {items.name}';
            alg   = [names(len1(1)+1:len1(2)),parameters(len1(1)+1:len1(2))];
            pro   = [names(len1(2)+1:len1(3)),parameters(len1(2)+1:len1(3))];
            % Decompose the parameter settings of problems
            pro1 = {};
            for i = 1 : size(pro,1)
                len = cellfun('length',pro{i,2});
                if any(len~=0&len~=1&len~=max(len))
                    errordlg(sprintf('The number of the parameters for <%s> is illegal (should be 0, 1 or %d).',pro{i,1},max(len)),get(gcf,'Name'),'modal'); beep;
                    return;
                else
                    subpro      = cell(max(1,max(len)),2);
                    subpro(:,1) = pro(i,1);
                    subpro(:,2) = {cell(1,length(pro{i,2}))};
                    for j = 1 : size(subpro,1)
                        subpro{j,2}(len>1) = num2cell(cellfun(@(S)S(j),pro{i,2}(len>1)));
                        subpro{j,2}(len<2) = pro{i,2}(len<2);
                    end
                    pro1 = [pro1;subpro];
                end
            end
            % Generate the setting of number of results
            pops = str2num(obj.control.numpopEdit.handle.String);
            if ~isa(pops,'double') || ~isreal(pops) || ~isscalar(pops) || pops~=fix(pops) || pops < 1
                errordlg('The number of results is illegal.',get(gcf,'Name'),'modal'); beep;
                return;
            end
            % Generate the setting of number of runs
            runtimes = str2num(obj.control.runtimeEdit.handle.String);
            if ~isa(runtimes,'double') || ~isreal(runtimes) || ~isscalar(runtimes) || runtimes~=fix(runtimes) || runtimes < 1
                errordlg('The number of runs is illegal.',get(gcf,'Name'),'modal'); beep;
                return;
            end
            % Save the current experimental settings
            if ~obj.cb_slSetting(obj.control.folderEdit.string)
                return;
            end
            % Initialise the data
            try
                obj.data = module_experiment_result(alg,pro1,runtimes,pops,fileparts(obj.control.folderEdit.handle.String));
            catch err
                errordlg(err.message,get(gcf,'Name'),'modal'); beep;
                return;
            end
            % Initialise the table
            obj.initialiseTable();
            % Update the state of controls
            obj.GUI.menu.enable           = false;
            obj.panel.panels(1).enable    = false;
            obj.panel.panels(2).enable    = false;
            obj.control.playStart.visible = false;
            obj.control.playPause.visible = true;
            obj.control.playStop.value    = false;
            % Initialise the menu of Wilcoxon rank sum test
            obj.control.ranksumMenu.del(1:length(obj.control.ranksumMenu.items));
            for i = 1 : size(obj.data.Algorithms,1)
                obj.control.ranksumMenu.add([],obj.data.Algorithms{i,1},'choosed',true,'callback',@obj.cb_ranksumMenu);
            end
            obj.control.ranksumMenu.items(end).value = true;
            if ~inparallel
                % Start the experiment
                for i = 1 : obj.data.nRuns
                    for j = 1 : size(obj.data.Problems,1)
                        for k = 1 : size(obj.data.Algorithms,1)
                            % Initialize the progress bar
                            obj.control.playbar.value = 0;
                            % Highlight the current grid
                            obj.control.table.selected(j,size(obj.control.table.handle.Data,2)-size(obj.data.Result,2)+k);
                            drawnow();
                            filename = obj.data.filename(k,j,i);
                            if exist(filename,'file') ~= 2
                                % Execute the algorithm
                                Inputs = {'-N',obj.data.Problems{j,2}{1},'-M',obj.data.Problems{j,2}{2},'-D',obj.data.Problems{j,2}{3},'-evaluation',obj.data.Problems{j,2}{4},...
                                          '-algorithm',[{str2func(obj.data.Algorithms{k,1})},obj.data.Algorithms{k,2}],...
                                          '-problem',[{str2func(obj.data.Problems{j,1})},obj.data.Problems{j,2}(5:end)],'-save',obj.data.nPops,'-outputFcn',@obj.outputFcn};
                                try
                                    Global = GLOBAL(Inputs{:},'-run',i);
                                    Global.Start();
                                    if ~obj.control.playStop.value
                                        obj.data.Result{j,k,i} = load(filename,'result','metric','-mat');
                                        obj.flushTable(j);
                                    end
                                catch err
                                    obj.cb_stop();
                                    rethrow(err);
                                end
                            else
                                % Load existing result
                                obj.data.Result{j,k,i} = load(filename,'result','metric','-mat');
                                obj.flushTable(j);
                            end
                            if obj.control.playStop.value
                                return;
                            end
                        end
                    end
                end
            else
                % Start the experiment in parallel
                for j = 1 : size(obj.data.Problems,1)
                    for k = 1 : size(obj.data.Algorithms,1)
                        % Initialize the progress bar
                        obj.control.playbar.value = 0;
                        obj.control.playlabel.handle.String = sprintf('<Parallel> %s_%s_M%d_D%d (%5.1f%%)',...
                                                              obj.data.Algorithms{k,1},obj.data.Problems{j,1},obj.data.Problems{j,3}{2},obj.data.Problems{j,3}{3},0);
                        % Highlight the current grid
                        obj.control.table.selected(j,size(obj.control.table.handle.Data,2)-size(obj.data.Result,2)+k);
                        drawnow();
                        % Execute the algorithm in parallel
                        Folder   = obj.data.Folder;
                        Filename = arrayfun(@(i)obj.data.filename(k,j,i),1:obj.data.nRuns,'UniformOutput',false);
                        Inputs   = {'-N',obj.data.Problems{j,2}{1},'-M',obj.data.Problems{j,2}{2},'-D',obj.data.Problems{j,2}{3},'-evaluation',obj.data.Problems{j,2}{4},...
                                    '-algorithm',[{str2func(obj.data.Algorithms{k,1})},obj.data.Algorithms{k,2}],...
                                    '-problem',[{str2func(obj.data.Problems{j,1})},obj.data.Problems{j,2}(5:end)],'-save',obj.data.nPops,'-outputFcn',@(g)outputFcn2(Folder,g)};
                        try
                            for i = 1 : obj.data.nRuns
                                F(i) = parfeval(@parRun,0,Filename{i},Inputs{:},'-run',i);
                            end
                            numFinish = 0;
                            while numFinish < obj.data.nRuns
                                drawnow();
                                if obj.control.playStop.value
                                    cancel(F);
                                    return;
                                end
                                if obj.control.playStart.visible
                                    waitfor(obj.control.playStart.handle,'UserData');
                                end
                                if obj.control.playStop.value
                                    cancel(F);
                                    return;
                                end
                                index = fetchNext(F,0.01);
                                if ~isempty(index)
                                    numFinish = numFinish + 1;
                                    obj.data.Result{j,k,index} = load(Filename{index},'result','metric','-mat');
                                    obj.flushTable(j);
                                    obj.control.playbar.value = numFinish./obj.data.nRuns;
                                    obj.control.playlabel.handle.String = sprintf('<Parallel> %s_%s_M%d_D%d (%5.1f%%)',...
                                                                          obj.data.Algorithms{k,1},obj.data.Problems{j,1},obj.data.Problems{j,3}{2},obj.data.Problems{j,3}{3},100*numFinish./obj.data.nRuns);
                                end
                            end
                        catch err
                            cancel(F);
                            obj.cb_stop();
                            rethrow(err);
                        end
                    end
                end
            end
            % Update the state of controls
            obj.cb_stop();
        end
        %% Change the basic algorithm for rank sum test
        function cb_ranksumMenu(obj,hObject,eventdata)
            [obj.control.ranksumMenu.items.value] = deal(false);
            hObject.value = true;
            obj.control.tableButton(7).value = true;
            obj.flushTable;
        end
        %% Right click on the table
        function cb_table(obj,hObject,eventdata)
            [row,column] = obj.control.table.selected;
            if row > 0 && row <= size(obj.data.Result,1) && column > size(obj.control.table.handle.Data,2)-size(obj.data.Result,2)
                obj.control.tableMenu.show();
            end
        end
        %% Show the figure of specified cell
        function cb_tableMenu(obj,type)
            [p,a]     = obj.control.table.selected;
            a         = a - size(obj.control.table.handle.Data,2) + size(obj.data.Result,2);
            metric    = obj.control.showButton.handle.String;
            minMetric = ~ismember(metric,obj.GUI.metList{1,2});
            existing  = find(~squeeze(cellfun('isempty',obj.data.Result(p,a,:))));
            if isempty(existing)
                return
            end
            if type == 5
                switch metric
                    case {'Run times','runtime'}
                        return;
                    otherwise
                        data(:,1) = [obj.data.Result{p,a,existing(1)}.result{:,1}]';
                        data(:,2) = mean(obj.data.metricValue(p,a,@obj.Metric,metric,1),1)';
                end
                figure; Draw(data,'-k.','LineWidth',1.5,'MarkerSize',10);
                title(sprintf('%s on %s',obj.data.Algorithms{a},obj.data.Problems{p}));
                xlabel('Number of evaluations');
                ylabel(metric);
            else
                switch metric
                    case 'Run times'
                        index = 1;
                    otherwise
                        [~,rank] = sort(obj.data.metricValue(p,a,@obj.Metric,metric,0));
                        if type == 1 || type == 3
                            index = rank(ceil(length(rank)/2));
                        elseif minMetric
                            index = rank(1);
                        else
                            index = rank(end);
                        end
                end
                if type == 1 || type == 2
                    data = obj.data.Result{p,a,(existing(index))}.result{end}.objs;
                else
                    data = obj.data.Result{p,a,(existing(index))}.result{end}.decs;
                end
                figure; Draw(data);
                title(sprintf('%s on %s',obj.data.Algorithms{a},obj.data.Problems{p}));
            end
        end
        %% Running - start
        function cb_start(obj,hObject,eventdata)
            obj.control.playPause.moved = true;
            if obj.panel.panels(1).enable
                obj.cb_run(0);
            else
                obj.control.playStart.visible         = false;
                obj.control.playPause.visible         = true;
                obj.control.playStart.handle.UserData = rand;
            end
        end
        %% Running - pause
        function cb_pause(obj,hObject,eventdata)
            obj.control.playStart.moved   = true;
            obj.control.playStart.visible = true;
            obj.control.playPause.visible = false;
        end
        %% Running - stop
        function cb_stop(obj,hObject,eventdata)
            if ~obj.panel.panels(1).enable
                obj.control.playStart.handle.UserData = rand;
                obj.GUI.menu.enable           = true;
                obj.panel.panels(1).enable    = true;
                obj.panel.panels(2).enable    = true;
                obj.control.playStart.visible = true;
                obj.control.playPause.visible = false;
            end
        end
        %% Initialise the table
        function initialiseTable(obj)
            RowName = obj.data.Problems(:,1);
            for i = length(RowName) : -1 : 2
                if strcmp(RowName{i},RowName{i-1})
                    RowName{i} = '';
                end
            end
            obj.control.table.handle.RowName = RowName;
            obj.control.table.handle.Data    = cell(size(obj.data.Problems,1),size(obj.data.Algorithms,1));
            obj.control.table.handle.Data(:) = {'-'};
            obj.updataTableColumn();
        end
        %% Show the specified columns
        function updataTableColumn(obj)
            if isempty(obj.data)
                return;
            end
            nP   = size(obj.data.Result,1);
            nA   = size(obj.data.Result,2);
            str  = get([obj.control.tableButton(2:5).handle],'String');
            show = [obj.control.tableButton(2:5).value];
            obj.control.table.handle.ColumnName   = [str(show);obj.data.Algorithms(:,1)];
            obj.control.table.handle.ColumnWidth  = [repmat({33},1,sum(show)),repmat({130},1,nA)];
            obj.control.table.handle.ColumnFormat = repmat({'char'},1,sum(show)+nA);
            temp = cellfun(@(S)S(show),obj.data.Problems(:,3),'UniformOutput',false);
            temp = cellfun(@(S)int2str(S),cat(1,temp{:}),'UniformOutput',false);
            obj.control.table.handle.Data = [[temp;repmat({''},size(obj.control.table.handle.Data,1)-nP,sum(show))],obj.control.table.handle.Data(:,end-nA+1:end)];
        end
        %% Flush the table
        function flushTable(obj,problemindex)
            if isempty(obj.data)
                return;
            end
            [nP,nA,nR] = size(obj.data.Result);
            % Lock the window
            if nargin < 2
                problemindex = 1 : nP;
                obj.GUI.figure.busy = true; % Can flush the figure
                obj.control.waitlabel.visible = true;
            end
            % Identify the metric
            obj.control.table.handle.Data(problemindex,end-nA+1:end) = {'-'};
            metric    = obj.control.showButton.handle.String;
            minMetric = ~ismember(metric,obj.GUI.metList{1,2});
            % Update the status of table menu
            [obj.control.tableMenu.items.enable] = deal(~any(strcmp(metric,{'Run times','runtime'})));
            % Calculate and show the metric values
            switch metric
                case 'Run times'
                    for p = problemindex
                        % Calculate the run times of each algorithm on each
                        % test instance
                        data  = sum(reshape(~cellfun('isempty',obj.data.Result(p,:,:)),nA,nR),2)';
                        valid = find(data>0);
                        for a = valid
                            obj.control.table.handle.Data{p,end-nA+a} = sprintf('%d',data(a));
                        end
                    end
                    obj.control.table.handle.RowName = obj.control.table.handle.RowName(1:nP);
                    obj.control.table.handle.Data    = obj.control.table.handle.Data(1:nP,:);
                otherwise
                    basic = find([obj.control.ranksumMenu.items.value],1);
                    for p = problemindex
                        % Ignore the data sets having no population
                        data  = sum(reshape(~cellfun('isempty',obj.data.Result(p,:,:)),nA,nR),2)';
                        valid = find(data>0);
                        % Calculate the metric values of all the results
                        % in this row
                        data      = cell(1,nA);
                        meanvalue = zeros(size(data));
                        stdvalue  = zeros(size(data));
                        for a = valid
                            if nargin < 2
                                drawnow();
                            end
                            data{a}      = obj.data.metricValue(p,a,@obj.Metric,metric,0);
                            meanvalue(a) = mean(data{a});
                            stdvalue(a)  = std(data{a});
                            if obj.control.tableButton(6).value
                                str = sprintf('%.4e (%.2e)',meanvalue(a),stdvalue(a));
                            else
                                str = sprintf('%.4e',meanvalue(a));
                            end
                            obj.control.table.handle.Data{p,end-nA+a} = strrep(strrep(str,'e-0','e-'),'e+0','e+');
                        end
                        % Ignore the data sets having metric value of NaN
                        valid(arrayfun(@(S)isnan(S),meanvalue(valid))) = [];
                        % Calculate the ranksum test results
                        if obj.control.tableButton(7).value && ~isempty(data{basic}) && ismember(basic,valid)
                            for a = valid
                                if a ~= basic
                                    [~,k] = ranksum(data{a},data{basic});
                                    if k
                                        if meanvalue(a)<meanvalue(basic)&&minMetric || meanvalue(a)>meanvalue(basic)&&~minMetric
                                            obj.control.table.handle.Data{p,end-nA+a} = [obj.control.table.handle.Data{p,end-nA+a},' +'];
                                        else
                                            obj.control.table.handle.Data{p,end-nA+a} = [obj.control.table.handle.Data{p,end-nA+a},' -'];
                                        end
                                    else
                                        obj.control.table.handle.Data{p,end-nA+a} = [obj.control.table.handle.Data{p,end-nA+a},' ='];
                                    end
                                end
                            end
                        end
                        % Detect the best result in each row
                        if ~isempty(valid)
                            if minMetric
                                [~,best] = min(meanvalue(valid));
                            else
                                [~,best] = max(meanvalue(valid));
                            end
                            obj.control.table.handle.Data{p,end-nA+valid(best)} = ['<html><Font color=#3333ff>',obj.control.table.handle.Data{p,end-nA+valid(best)}];
                        end
                    end
                    % Count the ranksum test results in each column
                    if obj.control.tableButton(7).value
                        if length(obj.control.table.handle.RowName) == nP
                            obj.control.table.handle.RowName = [obj.control.table.handle.RowName;'+/-/='];
                            obj.control.table.handle.Data    = [obj.control.table.handle.Data;repmat({''},1,size(obj.control.table.handle.Data,2))];
                        end
                        noEmpty  = cellfun(@(S)length(S)>1,obj.control.table.handle.Data(1:end-1,end-nA+1:end));
                        lastSign = cellfun(@(S)S(end),obj.control.table.handle.Data(1:end-1,end-nA+1:end));
                        nSign1   = sum(noEmpty & lastSign=='+',1);
                        nSign2   = sum(noEmpty & lastSign=='-',1);
                        nSign3   = sum(noEmpty & lastSign=='=',1);
                        for a = 1 : nA
                            if a ~= basic
                                obj.control.table.handle.Data{end,end-nA+a} = sprintf('%d/%d/%d',nSign1(a),nSign2(a),nSign3(a));
                            else
                                obj.control.table.handle.Data{end,end-nA+a} = '';
                            end
                        end
                    else
                        obj.control.table.handle.RowName = obj.control.table.handle.RowName(1:nP);
                        obj.control.table.handle.Data    = obj.control.table.handle.Data(1:nP,:);
                    end
            end
            % Unlock the window
            if nargin < 2
                obj.GUI.figure.busy = false;    % Can flush the figure
                obj.control.waitlabel.visible = false;
            end
        end
        %% Save the table
        function saveTable(obj)
            if ~isempty(obj.control.table.handle.Data)
                try
                    [Name,Path] = uiputfile({'*.xlsx';'*.tex'},'','new');
                    if ischar(Name)
                        [~,~,Type] = fileparts(Name);
                        obj.GUI.figure.busy = true;
                        switch Type
                            case '.xlsx'
                                obj.table2excel(fullfile(Path,Name),obj.control.showButton.handle.String,[{'Problem'},obj.control.table.handle.ColumnName';obj.control.table.handle.RowName,obj.control.table.handle.Data]);
                            case '.tex'
                                obj.table2latex(fullfile(Path,Name),[{'Problem'},obj.control.table.handle.ColumnName';obj.control.table.handle.RowName,obj.control.table.handle.Data]);
                        end
                        obj.GUI.figure.busy = false;
                    end
                catch err
                    obj.GUI.figure.busy = false;
                    errordlg('Fail to save the table due to unknown reasons.',get(gcf,'Name'),'modal');
                    rethrow(err);
                end
            end
        end
        %% Output function
        function outputFcn(obj,Global)
            obj.control.playbar.value = Global.evaluated/Global.evaluation;
            obj.control.playlabel.handle.String = sprintf('<Sequential> %s_%s_M%d_D%d_%d.mat (%5.1f%%)',...
                                                  func2str(Global.algorithm),class(Global.problem),Global.M,Global.D,Global.run,Global.evaluated/Global.evaluation*100);
            if obj.control.playStop.value
                error('GLOBAL:Termination','Algorithm has terminated');
            end
            if obj.control.playStart.visible
                waitfor(obj.control.playStart.handle,'UserData');
            end
            if obj.control.playStop.value
                error('GLOBAL:Termination','Algorithm has terminated');
            end
            if Global.evaluated >= Global.evaluation
                folder = fullfile(obj.data.Folder,func2str(Global.algorithm));
                [~,~]  = mkdir(folder);
                result         = Global.result;
                metric.runtime = Global.runtime;
                save(fullfile(folder,sprintf('%s_%s_M%d_D%d_%d.mat',func2str(Global.algorithm),class(Global.problem),Global.M,Global.D,Global.run)),'result','metric');
            end
        end
        %% Table to excel
        function table2excel(obj,filename,sheetname,Data)
            [x,y] = size(Data);
            % Convert the indices to Excel cell number
            function range = getRange(varargin)
                if nargin == 2
                    range = num2str(varargin{1});
                    while varargin{2} > 0
                        range = [char(65+mod(varargin{2}-1,26)),range];
                        varargin{2} = floor((varargin{2}-1)/26);
                    end
                else
                    range = [getRange(varargin{1:2}),':',getRange(varargin{3:4})];
                end
            end
            % Open the file and get the sheet
            try
                Excel = actxGetRunningServer('Excel.Application');
            catch
                Excel = actxserver('Excel.Application');
            end
            if exist(filename,'file')
                delete(filename);
            end
            Workbook = invoke(Excel.Workbooks,'Add');
            Workbook.SaveAs(filename);
            Sheet = Workbook.ActiveSheet;
            Sheet.Name = sheetname;
            % Set the column width
            head = y - size(obj.data.Result,2);
            Sheet.Range(getRange(1,1,x,1)).ColumnWidth = 10;
            if head >= 2
                Sheet.Range(getRange(1,2,x,head)).ColumnWidth = 6;
            end
            Sheet.Range(getRange(1,head+1,x,y)).ColumnWidth = 22;
            % Initialize the cells
            Range = Sheet.Range(getRange(1,1,x,y));
            Range.HorizontalAlignment = 3;
            Range.Font.Name = 'Times New Roman';
            if strcmp(Data{x,1},'+/-/=')
                Sheet.Range(getRange(x,1,x,head)).Merge;
                Sheet.Range(getRange(x,1,x,y)).NumberFormat = '@';
            end
            % Set the font color
            for i = 1 : x
                for j = 1 : y
                    if ~isempty(strfind(Data{i,j},'<html><Font color=#3333ff>'))
                        Sheet.Range(getRange(i,j)).Font.Color = 15282995;
                    end
                end
            end
            % Write the data
            Range.Value = regexprep(Data,'<html><Font color=#3333ff>','','ignorecase');
            % Set the border and merge the cells
            Range.Borders.LineStyle = 1;
            for i = 2 : x
                if isempty(Data{i,1})
                    Sheet.Rows.Item(i-1).Borders.Item(4).Linestyle = 0;
                    Sheet.Rows.Item(i).Borders.Item(3).Linestyle   = 0;
                    Sheet.Range(getRange(i-1,1,i,1)).Merge;
                end
            end
            % Close the file
            Workbook.Save;
            Workbook.Close;
            Excel.Quit;
            Excel.delete;
        end
        %% Table to LaTeX
        function table2latex(obj,filename,Data)
            nP = size(obj.data.Result,1);
            nA = size(obj.data.Result,2);
            % Convert the data
            mainData = Data(2:nP+1,end-nA+1:end);
            mainData = regexprep(mainData,'<html><Font color=#3333ff>','\\hl{','ignorecase');
            mainData = regexprep(mainData,'+$','$+$');
            mainData = regexprep(mainData,'-$','$-$');
            mainData = regexprep(mainData,'=$','$\\approx$');
            temp     = ~cellfun('isempty',strfind(mainData,'\hl{'));
            mainData(temp)            = strcat(mainData(temp),'}');
            Data(2:nP+1,end-nA+1:end) = mainData;
            Data(end,1)          = regexprep(Data(end,1),'^\+/\-/=$',['\\multicolumn{',num2str(size(Data,2)-nA),'}{c}{$+/-/\\approx$}']);
            Data(1,2:end-nA)     = strcat('$',Data(1,2:end-nA),'$');
            Data(1,end-nA+1:end) = regexprep(Data(1,end-nA+1:end),'_','\\_');
            noEmpty = ~cellfun('isempty',Data(:,1));
            for i = 2 : nP+1
                if noEmpty(i)
                    Data{i,1} = sprintf('\\multirow{%d}{*}{%s}',find([noEmpty(i+1:end);true],1),Data{i,1});
                end
            end
            % Generate the LaTeX code
            Code = eval(sprintf('strcat(%s)',strjoin(arrayfun(@(S)num2str(S,'Data(:,%d)'),1:size(Data,2),'UniformOutput',false),',''&'',')));
            Code = strcat(Code,'\\');
            if ~isempty(regexp(Code{end,1},'^\\multicolumn','once'))
                temp = strfind(Code{end,1},'&');
                Code{end,1}(temp(1:size(Data,2)-nA-1)) = [];
            end
            noEmpty = find(noEmpty);
            for i = 3 : length(noEmpty)
                Code = [Code(1:noEmpty(i)+i-4);'\hline';Code(noEmpty(i)+i-3:end)];
            end
            Code = ['\documentclass[journal]{IEEEtran}'
                    '\usepackage{multirow,booktabs,color,soul,threeparttable}'
                    '\definecolor{hl}{rgb}{0.75,0.75,0.75}'
                    '\sethlcolor{hl}'
                    '\begin{document}'
                    '\begin{table*}[htbp]'
                    '\renewcommand{\arraystretch}{1.2}'
                    '\centering'
                    '\caption{No Title}'
                    ['\begin{tabular}{',repmat('c',1,size(Data,2)),'}']
                    '\toprule'
                    Code(1)
                    '\midrule'
                    Code(2:end)
                    '\bottomrule'
                    '\end{tabular}'
                    '\label{No Label}'
                    '\end{table*}'
                    '\end{document}'];
            fid = fopen(filename,'wt');
            for i = 1 : length(Code)
                fprintf(fid,'%s\n',Code{i});
            end
            fclose(fid);
        end
    end
end

%% Function for parallel execution
function parRun(filename,varargin)
    if exist(filename,'file') ~= 2
        Global = GLOBAL(varargin{:});
        Global.Start();
    end
end

%% Output function for parallelization
function outputFcn2(folder,Global)
    if Global.evaluated >= Global.evaluation
        folder = fullfile(folder,func2str(Global.algorithm));
        [~,~]  = mkdir(folder);
        result         = Global.result;
        metric.runtime = Global.runtime;
        save(fullfile(folder,sprintf('%s_%s_M%d_D%d_%d.mat',func2str(Global.algorithm),class(Global.problem),Global.M,Global.D,Global.run)),'result','metric');
    end
end