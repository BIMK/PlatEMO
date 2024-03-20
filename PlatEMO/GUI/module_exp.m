classdef module_exp < handle
%module_exp - Experimental module.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        GUI;                % The GUI object
        app  = struct();	% All the components
        data = [];          % All the results
    end
    methods(Access = ?GUI)
        %% Constructor
        function obj = module_exp(GUI)
            % The main grid
            obj.GUI = GUI;
            obj.app.maingrid = GUI.APP(3,1,uigridlayout(obj.GUI.app.maingrid,'RowHeight',{20,'1x'},'ColumnWidth',{'1.2x',1,'1x',1,'4x'},'Padding',[0 5 0 5],'RowSpacing',5,'ColumnSpacing',0,'BackgroundColor','w'));
            obj.app.label(1) = GUI.APP(1,1,uilabel(obj.app.maingrid,'Text','Algorithm selection','HorizontalAlignment','center','FontSize',11,'FontColor',[.5 .5 .5]));
            obj.app.label(2) = GUI.APP(1,3,uilabel(obj.app.maingrid,'Text','Parameter setting','HorizontalAlignment','center','FontSize',11,'FontColor',[.5 .5 .5]));
            obj.app.label(3) = GUI.APP(1,5,uilabel(obj.app.maingrid,'Text','Result display','HorizontalAlignment','center','FontSize',11,'FontColor',[.5 .5 .5]));
            GUI.APP([1 2],2,uipanel(obj.app.maingrid,'BackgroundColor',[.8 .8 .8]));
            GUI.APP([1 2],4,uipanel(obj.app.maingrid,'BackgroundColor',[.8 .8 .8]));

            % The first panel
            obj.app.grid(1)    = GUI.APP(2,1,uigridlayout(obj.app.maingrid,'RowHeight',{16,19,16,19,19,16,19,19,19,20,'1x',20,'1x',21,21,21},'ColumnWidth',{'1x','1x','1x'},'Padding',[8 5 8 0],'RowSpacing',3,'ColumnSpacing',5,'BackgroundColor','w'));
            obj.app.labelA(1)  = GUI.APP(1,[1 3],uilabel(obj.app.grid(1),'Text','Number of objectives','VerticalAlignment','bottom','FontSize',12,'FontColor',[.15 .6 .2],'FontWeight','bold'));
            obj.app.stateA(1)  = GUI.APP(2,1,uibutton(obj.app.grid(1),'state','Text','single','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',0,'Tooltip','The problem has a single objective','ValueChangedFcn',{@obj.cb_filter,1}));
            obj.app.stateA(2)  = GUI.APP(2,2,uibutton(obj.app.grid(1),'state','Text','multi','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',1,'Tooltip','The problem has 2 or 3 objectives','ValueChangedFcn',{@obj.cb_filter,2}));
            obj.app.stateA(3)  = GUI.APP(2,3,uibutton(obj.app.grid(1),'state','Text','many','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',0,'Tooltip','The problem has more than 3 objectives','ValueChangedFcn',{@obj.cb_filter,3}));
            obj.app.labelA(2)  = GUI.APP(3,[1 3],uilabel(obj.app.grid(1),'Text','Encoding scheme','VerticalAlignment','bottom','FontSize',12,'FontColor',[.15 .6 .2],'FontWeight','bold'));
            obj.app.stateA(4)  = GUI.APP(4,1,uibutton(obj.app.grid(1),'state','Text','real','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',1,'Tooltip','The decision variables are real numbers','ValueChangedFcn',{@obj.cb_filter,4}));
            obj.app.stateA(5)  = GUI.APP(4,2,uibutton(obj.app.grid(1),'state','Text','integer','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',0,'Tooltip','The decision variables are integers','ValueChangedFcn',{@obj.cb_filter,5}));
            obj.app.stateA(6)  = GUI.APP(4,3,uibutton(obj.app.grid(1),'state','Text','label','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',0,'Tooltip','The decision variables are labels','ValueChangedFcn',{@obj.cb_filter,6}));
            obj.app.stateA(7)  = GUI.APP(5,1,uibutton(obj.app.grid(1),'state','Text','binary','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',0,'Tooltip','The decision variables are binary numbers','ValueChangedFcn',{@obj.cb_filter,7}));
            obj.app.stateA(8)  = GUI.APP(5,2,uibutton(obj.app.grid(1),'state','Text','permutation','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',0,'Tooltip','The decision vector is a permutation','ValueChangedFcn',{@obj.cb_filter,8}));
            obj.app.labelA(3)  = GUI.APP(6,[1 3],uilabel(obj.app.grid(1),'Text','Special difficulties','VerticalAlignment','bottom','FontSize',12,'FontColor',[.15 .6 .2],'FontWeight','bold'));
            obj.app.stateA(9)  = GUI.APP(7,1,uibutton(obj.app.grid(1),'state','Text','large','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',0,'Tooltip','The problem has more than 100 decision variables','ValueChangedFcn',{@obj.cb_filter,9}));
            obj.app.stateA(10) = GUI.APP(7,2,uibutton(obj.app.grid(1),'state','Text','constrained','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',0,'Tooltip','The problem has constraints','ValueChangedFcn',{@obj.cb_filter,10}));
            obj.app.stateA(11) = GUI.APP(7,3,uibutton(obj.app.grid(1),'state','Text','expensive','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',0,'Tooltip','The objectives are computationally time-consuming','ValueChangedFcn',{@obj.cb_filter,11}));
            obj.app.stateA(12) = GUI.APP(8,1,uibutton(obj.app.grid(1),'state','Text','multimodal','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',0,'Tooltip','The objectives are multimodal','ValueChangedFcn',{@obj.cb_filter,12}));
            obj.app.stateA(13) = GUI.APP(8,2,uibutton(obj.app.grid(1),'state','Text','sparse','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',0,'Tooltip','Most decision variables of the optimal solutions are zero','ValueChangedFcn',{@obj.cb_filter,13}));
            obj.app.stateA(14) = GUI.APP(8,3,uibutton(obj.app.grid(1),'state','Text','dynamic','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',0,'Tooltip','The objectives vary periodically','ValueChangedFcn',{@obj.cb_filter,14}));
            obj.app.stateA(15) = GUI.APP(9,1,uibutton(obj.app.grid(1),'state','Text','multitask','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',0,'Tooltip','The problem has multiple tasks to be solved simultaneously','ValueChangedFcn',{@obj.cb_filter,15}));
            obj.app.stateA(16) = GUI.APP(9,2,uibutton(obj.app.grid(1),'state','Text','bilevel','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',0,'Tooltip','The problem has two nested objectives','ValueChangedFcn',{@obj.cb_filter,16}));
            obj.app.stateA(17) = GUI.APP(9,3,uibutton(obj.app.grid(1),'state','Text','robust','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',0,'Tooltip','The objectives are influenced by uncertain factors','ValueChangedFcn',{@obj.cb_filter,17}));
            obj.app.labelA(4)  = GUI.APP(10,[1 2],uilabel(obj.app.grid(1),'Text','Algorithms','VerticalAlignment','bottom','FontSize',13,'FontColor',[.2 .4 .7],'FontWeight','bold'));
            obj.app.labelA(5)  = GUI.APP(10,3,uilabel(obj.app.grid(1),'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',10,'FontColor',[.2 .4 .7]));
            obj.app.listA(1)   = GUI.APP(11,[1 3],uilistbox(obj.app.grid(1),'FontColor',[.2 .4 .7],'ValueChangedFcn',{@obj.cb_updateList,2}));
            obj.app.labelA(6)  = GUI.APP(12,[1 2],uilabel(obj.app.grid(1),'Text','Problems','VerticalAlignment','bottom','FontSize',13,'FontColor',[.9 .5 .2],'FontWeight','bold'));
            obj.app.labelA(7)  = GUI.APP(12,3,uilabel(obj.app.grid(1),'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',10,'FontColor',[.9 .5 .2]));
            obj.app.listA(2)   = GUI.APP(13,[1 3],uilistbox(obj.app.grid(1),'FontColor',[.9 .5 .2],'ValueChangedFcn',{@obj.cb_updateList,-2}));
            obj.app.labelA(8)  = GUI.APP(14,[1 2],uilabel(obj.app.grid(1),'Text','Number of runs','FontColor',[.15 .6 .2],'FontWeight','bold','Tooltip','Number of runs for each algorithm on each problem'));
            obj.app.editA(1)   = GUI.APP(14,3,uieditfield(obj.app.grid(1),'numeric','Value',30,'limits',[1 inf],'RoundFractionalValues','on','Tooltip','Number of runs for each algorithm on each problem'));
            obj.app.labelA(9)  = GUI.APP(15,[1 2],uilabel(obj.app.grid(1),'Text','Number of results','FontColor',[.15 .6 .2],'FontWeight','bold','Tooltip','Number of populations saved in each run'));
            obj.app.editA(2)   = GUI.APP(15,3,uieditfield(obj.app.grid(1),'numeric','Value',1,'limits',[1 inf],'RoundFractionalValues','on','Tooltip','Number of populations saved in each run'));
            tempGrid           = GUI.APP(16,[1 3],uigridlayout(obj.app.grid(1),'RowHeight',{'1x'},'ColumnWidth',{'0.5x',20,'1x'},'Padding',[0 0 0 0],'RowSpacing',0,'ColumnSpacing',5,'BackgroundColor','w'));
            obj.app.labelA(10) = GUI.APP(1,1,uilabel(tempGrid,'Text','File path','FontColor',[.15 .6 .2],'FontWeight','bold','Tooltip','File path for saving experimental settings'));
            obj.app.buttonA    = GUI.APP(1,2,uibutton(tempGrid,'Text','...','BackgroundColor','w','ButtonpushedFcn',@obj.cb_filepath));
            obj.app.editA(3)   = GUI.APP(1,3,uieditfield(tempGrid,'Value',fullfile(cd,'Data','Setting.mat'),'HorizontalAlignment','right','Tooltip','File path for saving experimental settings'));
            
            % The second panel
            obj.app.listB   = uilist(obj.app.maingrid,obj.GUI.app.figure,obj.GUI.icon);
            obj.app.grid(2) = GUI.APP(2,3,obj.app.listB.grid);

            % The third panel
            obj.app.grid(3)    = GUI.APP(2,5,uigridlayout(obj.app.maingrid,'RowHeight',{30,'1x',30},'ColumnWidth',{230,'1x','1x',230},'Padding',[20 10 20 0],'ColumnSpacing',20,'BackgroundColor','w'));
            tempGrid           = GUI.APP(1,[1 4],uigridlayout(obj.app.grid(3),'RowHeight',{1,'1x',1},'ColumnWidth',{18,24,24,24,24,24,'1x','1x','1x','1.2x'},'Padding',[5 5 5 5],'RowSpacing',0,'ColumnSpacing',8,'BackgroundColor',[.95 .95 1]));
            tempPanel          = GUI.APP(2,1,uipanel(tempGrid,'BorderType','none','BackgroundColor',[.95 .95 1]));
            obj.app.toolC(1)   = uibutton(tempPanel,'Position',[-2.5 -2.5 24 24],'Text','','Icon',obj.GUI.icon.savetable,'BackgroundColor',[.95 .95 1],'Tooltip','Save the table','ButtonpushedFcn',@obj.cb_save);
            tempPanel          = GUI.APP(2,2,uipanel(tempGrid,'BorderType','none','BackgroundColor',[.95 .95 1]));
            obj.app.toolC(2)   = uibutton(tempPanel,'Position',[-2.5 -2.5 31 24],'Text','','Icon',obj.GUI.icon.figure,'BackgroundColor',[.95 .95 1],'Tooltip','Display the results','ButtonpushedFcn',@obj.cb_tableDisplay);
            obj.app.toolC(3)   = GUI.APP([1 3],3,uibutton(tempGrid,'state','Text','N','BackgroundColor',[.95 .95 1],'Tooltip','Show the population size','ValueChangedFcn',@obj.TableUpdateColumn));
            obj.app.toolC(4)   = GUI.APP([1 3],4,uibutton(tempGrid,'state','Text','M','BackgroundColor',[.95 .95 1],'Value',1,'Tooltip','Show the number of objectives','ValueChangedFcn',@obj.TableUpdateColumn));
            obj.app.toolC(5)   = GUI.APP([1 3],5,uibutton(tempGrid,'state','Text','D','BackgroundColor',[.95 .95 1],'Value',1,'Tooltip','Show the number of decision variables','ValueChangedFcn',@obj.TableUpdateColumn));
            obj.app.toolC(6)   = GUI.APP([1 3],6,uibutton(tempGrid,'state','Text','FE','BackgroundColor',[.95 .95 1],'Tooltip','Show the maximum number of function evaluations','ValueChangedFcn',@obj.TableUpdateColumn));
            obj.app.dropC(1)   = GUI.APP([1 3],7,uidropdown(tempGrid,'BackgroundColor',[.95 .95 1],'Tooltip','Show the specific metric values','Items',{},'Interruptible','off','BusyAction','cancel','ValueChangedFcn',@obj.TableUpdate));
            obj.app.dropC(2)   = GUI.APP([1 3],8,uidropdown(tempGrid,'BackgroundColor',[.95 .95 1],'Tooltip','Show the mean value and standard deviation','Items',{'Mean','Mean (STD)','Median','Median (IQR)'},'ItemsData',1:4,'Value',2,'Interruptible','off','BusyAction','cancel','ValueChangedFcn',@obj.TableUpdate));
            obj.app.dropC(3)   = GUI.APP([1 3],9,uidropdown(tempGrid,'BackgroundColor',[.95 .95 1],'Tooltip','Perform the Wilcoxon rank sum test','Items',{'none','Signed rank test','Rank sum test','Friedman test'},'ItemsData',1:4,'Value',3,'Interruptible','off','BusyAction','cancel','ValueChangedFcn',@obj.TableUpdate));
            obj.app.dropC(4)   = GUI.APP([1 3],10,uidropdown(tempGrid,'BackgroundColor',[.95 .95 1],'Tooltip','Highlight the best result','Items',{'none','Highlight the best','Highlight all the bests'},'ItemsData',1:3,'Value',2,'Interruptible','off','BusyAction','cancel','ValueChangedFcn',@obj.TableUpdate));
            obj.app.table      = GUI.APP(2,[1 4],uitable(obj.app.grid(3),'CellSelectionCallback',@obj.cb_tableSelect));
            obj.app.checkC     = GUI.APP(3,1,uicheckbox(obj.app.grid(3),'Text','Parallel execution','Tooltip','Perform the experiment with multiple CPUs','Value',0,'Enable',~isempty(ver('parallel'))));
            obj.app.buttonC(1) = GUI.APP(3,2,uibutton(obj.app.grid(3),'push','Text','Start','FontSize',16,'ButtonpushedFcn',@obj.cb_start));
            obj.app.buttonC(2) = GUI.APP(3,3,uibutton(obj.app.grid(3),'push','Text','Stop','FontSize',16,'Enable','off','ButtonpushedFcn',@obj.cb_stop));
            obj.app.labelC     = GUI.APP(3,4,uilabel(obj.app.grid(3),'Text','','HorizontalAlignment','right','VerticalAlignment','center'));
            obj.app.tablemenu  = uicontext(obj.GUI.app.figure,110);
            obj.app.tablemenu.add('  Populations (obj.)','',{@obj.cb_tableShow,1});
            obj.app.tablemenu.add('  Populations (dec.)','',{@obj.cb_tableShow,2});
            obj.app.tablemenu.add('  Metric values','',{@obj.cb_tableShow,3});
            obj.app.tablemenu.flush();
            
            % Initialization
            obj.cb_filter([],[],2);
        end
    end
    methods(Access = private)
        %% Update the algorithms and problems in the lists
        function cb_filter(obj,~,~,index)
            if index < 4
                [obj.app.stateA(1:3).Value] = deal(0);
                obj.app.stateA(index).Value = 1;
            elseif index < 9
                [obj.app.stateA(4:8).Value] = deal(0);
                obj.app.stateA(index).Value = 1;
            end
            filter = [obj.app.stateA.Value];
            func   = @(s)all(any(repmat([true,filter],size(s,1),1)&s,2)) && all((any(s(:,2:end),1)&filter)==filter);
            % Update the list of algorithms
            show   = cellfun(func,obj.GUI.algList(:,1));
            obj.app.listA(1).Items = ['(Open File)';obj.GUI.algList(show,2)];
            obj.app.listA(1).Value = {};
            obj.app.labelA(5).Text = sprintf('%d / %d',sum(show),length(show));
            % Update the list of problems
            show   = cellfun(func,obj.GUI.proList(:,1));
            obj.app.listA(2).Items = ['(Open File)';obj.GUI.proList(show,2)];
            obj.app.listA(2).Value = {};
            obj.app.labelA(7).Text = sprintf('%d / %d',sum(show),length(show));
            % Update the lists of metrics
            show   = cellfun(@(s)func(s(2:end,1:end-2)),obj.GUI.metList(:,1));
            obj.app.dropC(1).Items = ['Number of runs';'runtime';obj.GUI.metList(show,2)];
            obj.TableUpdate();
        end
        %% Update the parameter list
        function cb_updateList(obj,~,~,type)
            if type > 0
                filename = obj.app.listA(1).Value;
                tip = 'ALGORITHM';
            else
                filename = obj.app.listA(2).Value;
                tip = 'PROBLEM';
            end
            if contains(filename,'Open File')
                [file,path] = uigetfile({'*.m','MATLAB class'},'');
                if file ~= 0
                    try
                        filename = fullfile(path,file);
                        f   = fopen(filename);
                        str = fgetl(f);
                        fclose(f);
                        assert(contains(str,['< ',tip]));
                        addpath(path);
                    catch
                        uialert(obj.GUI.app.figure,sprintf('The selected file is not a subclass of %s.',tip),'Error');
                        return;
                    end
            	else
                    return;
                end
            else
                filename = [filename,'.m'];
            end
            obj.app.listB.add(filename,type);
            obj.app.listB.flush();
        end
        %% Load or save experimental settings
        function success = cb_filepath(obj,~,~,filename)
            success = false;
            if nargin < 4   % Load experimental settings
                [file,folder] = uigetfile({'*.mat','MAT file'},'',fileparts(obj.app.editA(3).Value));
                if ischar(file)
                    try
                        filename = fullfile(folder,file);
                        load(filename,'Setting','Environment','-mat');
                        obj.app.listB.del(1:length(obj.app.listB.items));
                        cellfun(@(s)obj.app.listB.add([s,'.m'],2),Setting{1});
                        cellfun(@(s)obj.app.listB.add([s,'.m'],-2),Setting{2});
                        set([obj.app.listB.items.edit],{'Value'},Setting{3});
                        obj.app.listB.flush();
                        set(obj.app.editA(1:2),{'Value'},Environment);
                        obj.app.editA(3).Value = filename;
                    catch
                        uialert(obj.GUI.app.figure,sprintf('Fail to load the experimental settings from %s.',filename),'Error');
                        return;
                    end
                end
            else            % Save experimental settings
                index       = find([obj.app.listB.items.type]<0,1);
                Setting{1}  = get([obj.app.listB.items(1:index-1).title],{'Text'});
                Setting{2}  = get([obj.app.listB.items(index:end).title],{'Text'});
                Setting{3}  = get([obj.app.listB.items.edit],{'Value'});
                Environment = get([obj.app.editA(1:2)],{'Value'});
                try
                    [folder,file] = fileparts(filename);
                    if isempty(file)
                        file = 'Setting';
                    end
                    filename = fullfile(folder,[file,'.mat']);
                    [~,~]    = mkdir(folder);
                    save(filename,'Setting','Environment','-mat');
                    obj.app.editA(3).Value = filename;
                catch
                    uialert(obj.GUI.app.figure,sprintf('Fail to save the experimental settings to %s.',filename),'Error');
                    return;
                end
            end
            success = true;
        end
        %% Start the execution
        function cb_start(obj,~,~)
            if strcmp(obj.app.buttonC(1).Text,'Pause')
                obj.app.buttonC(1).Text = 'Continue';
            elseif strcmp(obj.app.buttonC(1).Text,'Continue')
                obj.app.buttonC(1).Text = 'Pause';
            else
                try
                    % Check the validity of settings
                    isParallel = obj.app.checkC.Value;
                    ALG = [];
                    PRO = [];
                    assert(~isempty(obj.app.listB.items),'No algorithm is selected.');
                    type = [obj.app.listB.items.type];
                    assert(any(type>0),'No algorithm is selected.');
                    assert(any(type<0),'No problem is selected.');
                    allList  = [obj.GUI.algList;obj.GUI.proList];
                    allLabel = any(cell2mat(allList(ismember(allList(:,2),get([obj.app.listB.items.title],'Text')),1)),1);
                    SorM     = [allLabel(2),allLabel(3)||allLabel(4)];
                    assert(SorM(1)~=SorM(2),'Cannot perform single- and multi-objective optimization simultaneously.');
                    % Generate the ALGORITHM and PROBLEM object
                    for i = 1 : length(type)
                        item = obj.app.listB.items(i);
                        para = cell(1,length(item.edit));
                        for j = 1 : length(para)
                            if ~isempty(item.edit(j).Value)
                                para{j} = str2num(item.edit(j).Value);
                                assert(~isempty(para{j}),'The parameter "%s" of %s is illegal.',item.label(j).Text,item.title.Text);
                            end
                        end
                        if type(i) > 0
                            if ~isParallel
                                ALG = [ALG,feval(item.title.Text,'parameter',para,'save',obj.app.editA(2).Value,'outputFcn',@obj.outputFcn)];
                            else
                                ALG = [ALG,feval(item.title.Text,'parameter',para,'save',obj.app.editA(2).Value,'outputFcn',@(~,~)[])];
                            end
                        else
                            len = cellfun(@length,para(1:4));
                            for j = 1 : max(1,max(len))
                                paraSub        = para;
                                paraSub(len>1) = cellfun(@(s)s(min(end,j)),para(len>1),'UniformOutput',false);
                                PRO            = [PRO,feval(item.title.Text,'N',paraSub{1},'M',paraSub{2},'D',paraSub{3},item.label(4).Text,paraSub{4},'parameter',paraSub(5:end))];
                            end
                        end
                    end
                catch err
                    uialert(obj.GUI.app.figure,err.message,'Invalid parameter settings');
                    return;
                end
                % Save the experimental settings
                if ~obj.cb_filepath([],[],obj.app.editA(3).Value)
                    return;
                end
                obj.data = struct('ALG',ALG,'PRO',PRO,'folder',fileparts(obj.app.editA(3).Value),'result',{cell(length(PRO),length(ALG),obj.app.editA(1).Value)},'metric',{cell(length(PRO),length(ALG),obj.app.editA(1).Value)});
                % Initialize the table
                rowName = arrayfun(@class,PRO,'UniformOutput',false);
                for i = length(rowName) : -1 : 2
                    if strcmp(rowName{i},rowName{i-1})
                        rowName{i} = '';
                    end
                end
                obj.app.table.RowName = rowName;
                obj.app.table.Data    = cell(length(PRO),length(ALG));
                obj.TableUpdateColumn();
                % Update the GUI
                [obj.GUI.app.button.Enable] = deal('off');
                [obj.app.stateA.Enable]     = deal('off');
                [obj.app.listA.Enable]      = deal('off');
                [obj.app.editA(1:2).Enable] = deal('off');
                [obj.app.editA(3).Enable]   = deal('off');
                obj.app.buttonA.Enable      = deal('off');
                obj.app.listB.Enable        = deal('off');
                obj.app.buttonC(1).Text     = 'Pause';
                obj.app.buttonC(2).Enable   = 'on';
                drawnow('limitrate');
                % Perform the experiment
                for p = 1 : length(PRO)
                    for a = 1 : length(ALG)
                        % Load existing results
                        arrayfun(@(r)obj.ResultLoad(p,a,r),1:size(obj.data.result,3));
                        obj.TableUpdate([],[],p);
                        runIndex = find(reshape(cellfun(@isempty,obj.data.result(p,a,:)),1,[]));
                        % Run algorithms in sequence
                        if ~isempty(runIndex) && ~isParallel
                            try
                                for r = runIndex
                                    ALG(a).Solve(PRO(p));
                                    obj.ResultSave(p,a,r,ALG(a).result,ALG(a).metric);
                                    obj.ResultLoad(p,a,r);
                                    obj.TableUpdate([],[],p);
                                    if strcmp(obj.app.buttonC(2).Enable,'off')
                                        return;
                                    end
                                end
                            catch err
                                uialert(obj.GUI.app.figure,'The algorithm terminates unexpectedly, please refer to the command window for details.','Error');
                                obj.cb_stop();
                                rethrow(err);
                            end
                        % Run algorithms in parallel
                        elseif ~isempty(runIndex)
                            try
                                Future = arrayfun(@(s)parfeval(@parallelFcn,2,ALG(a),PRO(p)),runIndex);
                                while ~all([Future.Read])
                                    drawnow('limitrate');
                                    if strcmp(obj.app.buttonC(2).Enable,'off')
                                        cancel(Future);
                                        return;
                                    end
                                    if strcmp(obj.app.buttonC(1).Text,'Continue')
                                        waitfor(obj.app.buttonC(1),'Text');
                                    end
                                    if strcmp(obj.app.buttonC(2).Enable,'off')
                                        cancel(Future);
                                        return;
                                    end
                                    [r,result,metric] = fetchNext(Future,0.01);
                                    if ~isempty(r)
                                        obj.ResultSave(p,a,runIndex(r),result,metric);
                                        obj.ResultLoad(p,a,runIndex(r));
                                        obj.TableUpdate([],[],p);
                                    end
                                end
                            catch err
                                try
                                    cancel(Future);
                                catch
                                end
                                uialert(obj.GUI.app.figure,'The algorithm terminates unexpectedly, please refer to the command window for details.','Error');
                                obj.cb_stop();
                                rethrow(err);
                            end
                        end
                    end
                end
                obj.cb_stop();
            end
        end
        %% Stop the execution
        function cb_stop(obj,~,~)
            [obj.GUI.app.button.Enable] = deal('on');
            [obj.app.stateA.Enable]     = deal('on');
            [obj.app.listA.Enable]      = deal('on');
            [obj.app.editA(1:2).Enable] = deal('on');
            [obj.app.editA(3).Enable]   = deal('on');
            obj.app.buttonA.Enable      = deal('on');
            obj.app.listB.Enable        = deal('on');
            obj.app.buttonC(1).Text     = 'Start';
            obj.app.buttonC(2).Enable   = 'off';
        end
        %% Output function
        function outputFcn(obj,Algorithm,Problem)
            assert(strcmp(obj.app.buttonC(2).Enable,'on'),'PlatEMO:Termination','');
            if strcmp(obj.app.buttonC(1).Text,'Continue')
                waitfor(obj.app.buttonC(1),'Text');
            end
            assert(strcmp(obj.app.buttonC(2).Enable,'on'),'PlatEMO:Termination','');
        end
        %% Show the specified columns
        function TableUpdateColumn(obj,~,~)
            if ~isempty(obj.data)
                % Update the columns
                str  = get(obj.app.toolC(3:6),'Text');
                show = [obj.app.toolC(3:6).Value];
                obj.app.table.ColumnName   = [str(show)',arrayfun(@class,obj.data.ALG,'UniformOutput',false)];
                obj.app.table.ColumnWidth  = [repmat({45},1,sum(show)),repmat({'auto'},1,length(obj.data.ALG))];
                obj.app.table.ColumnFormat = repmat({'char'},1,sum(show)+length(obj.data.ALG));
                oldsize = size(obj.app.table.Data,2);
                str = arrayfun(@num2str,[obj.data.PRO.N;obj.data.PRO.M;obj.data.PRO.D;obj.data.PRO.maxFE]','UniformOutput',false);
                obj.app.table.Data = [[str(:,show);repmat({''},size(obj.app.table.Data,1)-length(obj.data.PRO),sum(show))],obj.app.table.Data(:,end-length(obj.data.ALG)+1:end)];
                % Update the styles
                styleLoc = cat(1,obj.app.table.StyleConfigurations.TargetIndex{:});
                if ~isempty(styleLoc)
                    obj.app.table.removeStyle();
                    obj.app.table.addStyle(uistyle('FontWeight','bold'),'cell',[styleLoc(:,1),styleLoc(:,2)+size(obj.app.table.Data,2)-oldsize]);
                end
            end
        end
        %% Update the table
        function TableUpdate(obj,~,~,proindex)
            % Change the tooltips of drop-down components
            str = {'Show the mean value','Show the mean value and standard deviation','Show the median value','Show the median value and interquartile range'};
            obj.app.dropC(2).Tooltip = str(obj.app.dropC(2).Value);
            str = {'','Perform the Wilcoxon signed rank test','Perform the Wilcoxon rank sum test','Perform the Friedman''s test'};
            obj.app.dropC(3).Tooltip = str(obj.app.dropC(3).Value);
            str = {'','Highlight the best result','Highlight the best and statistically similar results'};
            obj.app.dropC(4).Tooltip = str(obj.app.dropC(4).Value);
            if ~isempty(obj.data)
                [nP,nA,nR] = size(obj.data.result);
                if nargin < 4
                    proindex = 1 : nP;
                end
                % Delete the cells and styles in the table
                obj.app.table.Data(proindex,end-nA+1:end) = {''};
                styleLoc = cat(1,obj.app.table.StyleConfigurations.TargetIndex{:});
                if ~isempty(styleLoc)
                    styleLoc(ismember(styleLoc(:,1),proindex),:) = [];
                end
                obj.app.table.removeStyle();
                % Identify the metric
                metric = obj.app.dropC(1).Value;
                allMet = [obj.GUI.metList;{[0 1]},'Number of runs';{[1 0]},'runtime'];
                minMet = allMet{find(ismember(allMet(:,2),metric),1),1}(1,end-1);
                if strcmp(metric,'Number of runs')
                    for p = proindex
                        % Show the number of runs
                        cdata = sum(reshape(~cellfun(@isempty,obj.data.result(p,:,:)),nA,nR),2);
                        obj.app.table.Data(p,end-nA+1:end) = arrayfun(@num2str,cdata','UniformOutput',false);
                    end
                    % Hide the row of statistical results
                    obj.app.table.RowName = obj.app.table.RowName(1:nP);
                    obj.app.table.Data    = obj.app.table.Data(1:nP,:);
                else
                    for p = proindex
                        % Show the metric values
                        valid = find(any(reshape(~cellfun(@isempty,obj.data.result(p,:,:)),nA,nR),2))';
                        cdata = cell(1,nA);     % All metric values
                        mdata = zeros(1,nA);    % Mean or median value
                        sdata = zeros(1,nA);    % STD or IQR value
                        for a = valid
                            if nargin < 4
                                drawnow('limitrate');
                            end
                            cdata{a} = obj.GetMetricValue(p,a,metric,false);
                            datapure = cdata{a}(~isnan(cdata{a}));
                            if obj.app.dropC(2).Value < 3
                                mdata(a) = mean(datapure);
                                sdata(a) = std(datapure);
                            else
                                mdata(a) = median(datapure);
                                sdata(a) = iqr(datapure);
                            end
                            if ismember(obj.app.dropC(2).Value,[1 3])
                                str = sprintf('%.4e',mdata(a));
                            else
                                str = sprintf('%.4e (%.2e)',mdata(a),sdata(a));
                            end
                            obj.app.table.Data{p,end-nA+a} = strrep(strrep(str,'e-0','e-'),'e+0','e+');
                        end
                        % Highlight the best metric value
                        valid(arrayfun(@(s)isnan(s),mdata(valid))) = [];
                        if ~isempty(valid) && obj.app.dropC(4).Value > 1
                            if minMet
                                [~,best] = min(mdata(valid));
                            else
                                [~,best] = max(mdata(valid));
                            end
                            styleLoc = [styleLoc;p,size(obj.app.table.Data,2)-nA+valid(best)];
                        end
                        % Calculate the statistical test results
                        if obj.app.dropC(3).Value > 1 && length(valid) > 1 && ismember(nA,valid)
                            minlen = min(cellfun(@length,cdata(valid)));
                            if minlen > 1
                                vdata = cellfun(@(s)s(1:minlen),cdata(valid),'UniformOutput',false);
                                if obj.app.dropC(3).Value == 4
                                    [~,~,stats] = friedman(cell2mat(vdata),1,'off');
                                    c     = multcompare(stats,'Display','off');
                                    diff1 = c(any(c==length(vdata),2),end) < 0.05;
                                    if obj.app.dropC(4).Value == 3
                                        diff2 = c(any(c==best,2),end) < 0.05;
                                    end
                                elseif obj.app.dropC(3).Value == 3
                                    diff1 = cellfun(@(s)ranksum(s,vdata{end})<0.05,vdata(:,1:end-1));
                                    if obj.app.dropC(4).Value == 3
                                        diff2 = cellfun(@(s)ranksum(s,vdata{best})<0.05,vdata(:,[1:best-1,best+1:end]));
                                    end
                                elseif obj.app.dropC(3).Value == 2
                                    diff1 = cellfun(@(s)signrank(s,vdata{end})<0.05,vdata(:,1:end-1));
                                    if obj.app.dropC(4).Value == 3
                                        diff2 = cellfun(@(s)signrank(s,vdata{best})<0.05,vdata(:,[1:best-1,best+1:end]));
                                    end
                                end
                                for a = 1 : length(diff1)
                                    if ~diff1(a) || mdata(valid(a))==mdata(nA)
                                        obj.app.table.Data{p,end-nA+valid(a)} = [obj.app.table.Data{p,end-nA+valid(a)},' ='];
                                    elseif mdata(valid(a))<mdata(nA)&&minMet || mdata(valid(a))>mdata(nA)&&~minMet
                                        obj.app.table.Data{p,end-nA+valid(a)} = [obj.app.table.Data{p,end-nA+valid(a)},' +'];
                                    else
                                        obj.app.table.Data{p,end-nA+valid(a)} = [obj.app.table.Data{p,end-nA+valid(a)},' -'];
                                    end
                                    if obj.app.dropC(4).Value == 3 && ~diff2(a)
                                        styleLoc = [styleLoc;p,size(obj.app.table.Data,2)-nA+valid(a)+(a>=best)];
                                    end
                                end
                            end
                        end
                    end
                    if ~isempty(styleLoc)
                        obj.app.table.addStyle(uistyle('FontWeight','bold'),'cell',styleLoc);
                    end
                    % Count the statistical test results
                    if obj.app.dropC(3).Value > 1
                        if length(obj.app.table.RowName) == nP
                            obj.app.table.RowName = [obj.app.table.RowName;'+/-/='];
                            obj.app.table.Data    = [obj.app.table.Data;repmat({''},1,size(obj.app.table.Data,2))];
                        end
                        sign1 = cellfun(@(s)~isempty(s)&&strcmp('+',s(end)),obj.app.table.Data(1:end-1,end-nA+1:end));
                        sign2 = cellfun(@(s)~isempty(s)&&strcmp('-',s(end)),obj.app.table.Data(1:end-1,end-nA+1:end));
                        sign3 = cellfun(@(s)~isempty(s)&&strcmp('=',s(end)),obj.app.table.Data(1:end-1,end-nA+1:end));
                        for a = 1 : nA-1
                        	obj.app.table.Data{end,end-nA+a} = sprintf('%d/%d/%d',sum(sign1(:,a)),sum(sign2(:,a)),sum(sign3(:,a)));
                        end
                    else
                        obj.app.table.RowName = obj.app.table.RowName(1:nP);
                        obj.app.table.Data    = obj.app.table.Data(1:nP,:);
                    end
                end
            end
        end
        %% Load the result file
        function ResultLoad(obj,p,a,r)
            try
                filename = fullfile(obj.data.folder,class(obj.data.ALG(a)),sprintf('%s_%s_M%d_D%d_%d.mat',class(obj.data.ALG(a)),class(obj.data.PRO(p)),obj.data.PRO(p).M,obj.data.PRO(p).D,r));
                load(filename,'-mat','result','metric');
                obj.data.result{p,a,r} = result;
                obj.data.metric{p,a,r} = metric;
            catch
            end
        end
        %% Save the result file
        function ResultSave(obj,p,a,r,result,metric)
            folder   = fullfile(obj.data.folder,class(obj.data.ALG(a)));
            [~,~]    = mkdir(folder);
            filename = fullfile(folder,sprintf('%s_%s_M%d_D%d_%d.mat',class(obj.data.ALG(a)),class(obj.data.PRO(p)),obj.data.PRO(p).M,obj.data.PRO(p).D,r));
            save(filename,'result','metric');
        end
        %% Get the metric value 
        function score = GetMetricValue(obj,p,a,metName,showAll)
            metName = strrep(metName,' ','');
            score   = [];
            fes     = [];
            for r = find(reshape(~cellfun(@isempty,obj.data.result(p,a,:)),1,[]))
                if ~isfield(obj.data.metric{p,a,r},metName)
                    obj.data.metric{p,a,r}.(metName) = cellfun(@(S)obj.data.PRO(p).CalMetric(metName,S),obj.data.result{p,a,r}(:,2));
                    metric = obj.data.metric{p,a,r};
                    save(fullfile(obj.data.folder,class(obj.data.ALG(a)),sprintf('%s_%s_M%d_D%d_%d.mat',class(obj.data.ALG(a)),class(obj.data.PRO(p)),obj.data.PRO(p).M,obj.data.PRO(p).D,r)),'metric','-append');
                end
                try
                    if showAll  % Mean convergence profile of the metric
                        score = [score,obj.data.metric{p,a,r}.(metName)];
                        fes   = [fes,cell2mat(obj.data.result{p,a,r}(:,1))];
                    else        % All metric values of the last populations
                        score = [score;obj.data.metric{p,a,r}.(metName)(end)];
                    end
                catch
                end
            end
            if showAll
                score = [mean(fes,2),mean(score,2)];
            end
        end
        %% Save the table
        function cb_save(obj,~,~)
            if ~isempty(obj.app.table.Data)
                try
                    [Name,Path] = uiputfile({'*.xlsx','Excel table';'*.tex','TeX table';'*.txt','Text file';'*.mat','MAT file'},'','new');
                    if ischar(Name)
                        [~,~,Type] = fileparts(Name);
                        switch Type
                            case '.xlsx'
                                table2excel(fullfile(Path,Name),obj.app.dropC(1).Value,[{'Problem'},obj.app.table.ColumnName';obj.app.table.RowName,obj.app.table.Data],cat(1,obj.app.table.StyleConfigurations.TargetIndex{:}),size(obj.data.result,2));
                            case '.tex'
                                table2tex(fullfile(Path,Name),[{'Problem'},obj.app.table.ColumnName';obj.app.table.RowName,obj.app.table.Data],cat(1,obj.app.table.StyleConfigurations.TargetIndex{:}),size(obj.data.result,1),size(obj.data.result,2));
                            case '.txt'
                                table2txt(fullfile(Path,Name),[{'Problem'},obj.app.table.ColumnName';obj.app.table.RowName,obj.app.table.Data]);
                        	case '.mat'
                                table2mat(fullfile(Path,Name),obj.data.metric,obj.app.dropC(1).Value);
                        end
                    end
                catch err
                    uialert(obj.GUI.app.figure,'Fail to save the table, please refer to the command window for details.','Error');
                    rethrow(err);
                end
            end
        end
        %% Select the cells of table
        function cb_tableSelect(obj,~,event)
            if ~isempty(obj.data)
                grids = [event.Indices(:,1),event.Indices(:,2)-size(obj.app.table.Data,2)+size(obj.data.result,2)];
                grids(grids(:,1)>size(obj.data.result,1)|grids(:,2)<1,:) = [];
                if ~isempty(grids) && ~ismember(obj.app.dropC(1).Value,{'Number of runs','runtime'})
                    obj.app.table.UserData = [min(grids,[],1),max(grids,[],1)];
                else
                    obj.app.table.UserData = [];
                end
            end
        end
        %% Show the menu of result display
        function cb_tableDisplay(obj,~,~)
            if ~isempty(obj.app.table.UserData)
                obj.app.tablemenu.show();
            else
                uialert(obj.GUI.app.figure,'Please select at least one cell of metric value in the table.','Error');
            end
        end
        %% Show the results in new figure
        function cb_tableShow(obj,ui,~,type)
            ui.Parent.Visible = false;
            metric = obj.app.dropC(1).Value;
            loc    = obj.app.table.UserData;
            nRow   = loc(3) - loc(1) + 1;
            nCol   = loc(4) - loc(2) + 1;
            if type < 3
                movegui(figure('NumberTitle','off','Name','','Position',[0 0 240*nCol 220*nRow]),'center');
                for r = 1 : nRow
                    for c = 1 : nCol
                        ax    = Draw(axes('Unit','pixels','Position',[(c-1)*240+35 (nRow-r)*220+35 200 170]));
                        valid = find(reshape(~cellfun(@isempty,obj.data.result(r+loc(1)-1,c+loc(2)-1,:)),1,[]));
                        if ~isempty(valid)
                            [~,rank] = sort(obj.GetMetricValue(r+loc(1)-1,c+loc(2)-1,metric,false));
                            if type == 1
                                obj.data.PRO(r+loc(1)-1).DrawObj(obj.data.result{r+loc(1)-1,c+loc(2)-1,valid(rank(ceil(end/2)))}{end});
                            elseif type == 2
                                obj.data.PRO(r+loc(1)-1).DrawDec(obj.data.result{r+loc(1)-1,c+loc(2)-1,valid(rank(ceil(end/2)))}{end});
                            end
                        end
                        set(ax,'FontSize',8);
                        set(ax.Children,{'MarkerSize','LineWidth'},{4,0.6});
                        title([class(obj.data.ALG(c+loc(2)-1)),' on ',class(obj.data.PRO(r+loc(1)-1))],'Interpreter','none');
                    end
                end
            else
                movegui(figure('NumberTitle','off','Name','','Position',[0 0 290*nRow 260]),'center');
                for r = 1 : nRow
                    ax = Draw(axes('Unit','pixels','Position',[(r-1)*290+50 40 220 200]));
                    s  = {'o','+','s','*','^','x';'-k','--k','-b','--b','-g','--g'};
                    for c = 1 : nCol
                        value = obj.GetMetricValue(r+loc(1)-1,c+loc(2)-1,metric,true);
                        Draw(value,[s{1,mod(c-1,size(s,2))+1},s{2,mod(ceil(c/size(s,2))-1,size(s,2))+1}],'MarkerSize',5,'LineWidth',0.6);
                    end
                    legend(ax,arrayfun(@(s)class(s),obj.data.ALG(loc(2):loc(4)),'UniformOutput',false),'Location','best');
                    set(ax,'FontSize',10);
                    [ax.XLabel.String,ax.YLabel.String,ax.Title.String,ax.Title.Interpreter] = deal('Number of function evaluations',strrep(metric,'_',' '),class(obj.data.PRO(r+loc(1)-1)),'none');
                end
            end
        end
    end
end

%% Function for parallelization
function [result,metric] = parallelFcn(Algorithm,Problem)
    Algorithm.Solve(Problem);
    result = Algorithm.result;
    metric = Algorithm.metric;
end

%% Save the table to Excel
function table2excel(filename,sheetname,Data,styleLoc,nA)
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
    head = y - nA;
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
    for i = 1 : size(styleLoc,1)
        Sheet.Range(getRange(styleLoc(i,1)+1,styleLoc(i,2)+1)).Font.Color = 15282995;
    end
    % Write the data
    Range.Value = Data;
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

%% Save the table to TeX
function table2tex(filename,Data,styleLoc,nP,nA)
    % Convert the data
    mainData = Data(2:nP+1,end-nA+1:end);
    mainData = regexprep(mainData,'+$','$+$');
    mainData = regexprep(mainData,'-$','$-$');
    mainData = regexprep(mainData,'=$','$\\approx$');
    for i = 1 : size(styleLoc,1)
        mainData{styleLoc(i,1),styleLoc(i,2)-(size(Data,2)-nA-1)} = ['\hl{',mainData{styleLoc(i,1),styleLoc(i,2)-(size(Data,2)-nA-1)},'}'];
    end
    Data(2:nP+1,end-nA+1:end) = mainData;
    Data(end,1)          = regexprep(Data(end,1),'^\+/\-/=$',['\\multicolumn{',num2str(size(Data,2)-nA),'}{c}{$+/-/\\approx$}']);
    Data(1,2:end-nA)     = strcat('$',Data(1,2:end-nA),'$');
    Data(1,end-nA+1:end) = regexprep(Data(1,end-nA+1:end),'_','\\_');
    noEmpty = ~cellfun(@isempty,Data(:,1));
    for i = 2 : nP+1
        if noEmpty(i)
            Data{i,1} = sprintf('\\multirow{%d}{*}{%s}',find([noEmpty(i+1:end);true],1),Data{i,1});
        end
    end
    % Generate the TeX code
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

%% Save the table to .txt
function table2txt(filename,Data)
    fid = fopen(filename,'wt');
    for i = 1 : size(Data,1)
        fprintf(fid,'%s\n',strjoin(Data(i,:),'\t'));
    end
    fclose(fid);
end

%% Save the table to .mat
function table2mat(filename,Metric,metricName)
    metricName = strrep(metricName,' ','');
    if strcmp(metricName,'Numberofruns')
        Data = sum(~cellfun(@isempty,Metric),3);
    else
        Data = nan(size(Metric));
        for i = 1 : numel(Metric)
            if isfield(Metric{i},metricName)
                Data(i) = Metric{i}.(metricName)(end);
            end
        end
    end
    eval([metricName,'=Data;'])
    save(filename,metricName);
end