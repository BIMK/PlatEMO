classdef module_test < handle
%module_test - Test module.

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
        app  = struct(); 	% All the components
        data = {};          % All the results
    end
    methods(Access = ?GUI)
        %% Constructor
        function obj = module_test(GUI)
            % The main grid
            obj.GUI = GUI;
            obj.app.maingrid = GUI.APP(3,1,uigridlayout(obj.GUI.app.maingrid,'RowHeight',{20,'1x'},'ColumnWidth',{'1.2x',1,'1x',1,'2x','0.8x','0.2x',1,'1x'},'Padding',[0 5 0 5],'RowSpacing',5,'ColumnSpacing',0,'BackgroundColor','w'));
            obj.app.label(1) = GUI.APP(1,1,uilabel(obj.app.maingrid,'Text','Algorithm selection','HorizontalAlignment','center','FontSize',11,'FontColor',[.5 .5 .5]));
            obj.app.label(2) = GUI.APP(1,3,uilabel(obj.app.maingrid,'Text','Parameter setting','HorizontalAlignment','center','FontSize',11,'FontColor',[.5 .5 .5]));
            obj.app.label(3) = GUI.APP(1,[5 7],uilabel(obj.app.maingrid,'Text','Result display','HorizontalAlignment','center','FontSize',11,'FontColor',[.5 .5 .5]));
            obj.app.label(4) = GUI.APP(1,9,uilabel(obj.app.maingrid,'Text','Result selection','HorizontalAlignment','center','FontSize',11,'FontColor',[.5 .5 .5]));
            GUI.APP([1 2],2,uipanel(obj.app.maingrid,'BackgroundColor',[.8 .8 .8]));
            GUI.APP([1 2],4,uipanel(obj.app.maingrid,'BackgroundColor',[.8 .8 .8]));
            GUI.APP([1 2],8,uipanel(obj.app.maingrid,'BackgroundColor',[.8 .8 .8]));

            % The first panel
            obj.app.grid(1)    = GUI.APP(2,1,uigridlayout(obj.app.maingrid,'RowHeight',{16,19,16,19,19,16,19,19,19,22,'1x',22,'1x'},'ColumnWidth',{'1x','1x','1x'},'Padding',[8 10 8 0],'RowSpacing',3,'ColumnSpacing',5,'BackgroundColor','w'));
            obj.app.labelA(1)  = GUI.APP(1,[1 3],uilabel(obj.app.grid(1),'Text','Number of objectives','VerticalAlignment','bottom','FontSize',12,'FontColor',[.15 .6 .2],'FontWeight','bold'));
            obj.app.stateA(1)  = GUI.APP(2,1,uibutton(obj.app.grid(1),'state','Text','single','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',0,'Tooltip','The problem has a single objective','ValueChangedFcn',{@obj.cb_filter,1}));
            obj.app.stateA(2)  = GUI.APP(2,2,uibutton(obj.app.grid(1),'state','Text','multi','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',1,'Tooltip','The problem has 2 or 3 objectives','ValueChangedFcn',{@obj.cb_filter,2}));
            obj.app.stateA(3)  = GUI.APP(2,3,uibutton(obj.app.grid(1),'state','Text','many','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',0,'Tooltip','The problem has 4 or more objectives','ValueChangedFcn',{@obj.cb_filter,3}));
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
            obj.app.listA(1)   = GUI.APP(11,[1 3],uilistbox(obj.app.grid(1),'FontColor',[.2 .4 .7],'ValueChangedFcn',{@obj.cb_updateList,1}));
            obj.app.labelA(6)  = GUI.APP(12,[1 2],uilabel(obj.app.grid(1),'Text','Problems','VerticalAlignment','bottom','FontSize',13,'FontColor',[.9 .5 .2],'FontWeight','bold'));
            obj.app.labelA(7)  = GUI.APP(12,3,uilabel(obj.app.grid(1),'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',10,'FontColor',[.9 .5 .2]));
            obj.app.listA(2)   = GUI.APP(13,[1 3],uilistbox(obj.app.grid(1),'FontColor',[.9 .5 .2],'ValueChangedFcn',{@obj.cb_updateList,-1}));

            % The second panel
            obj.app.listB   = uilist(obj.app.maingrid,obj.GUI.app.figure,obj.GUI.icon);
            obj.app.grid(2) = GUI.APP(2,3,obj.app.listB.grid);

            % The third panel
            obj.app.dropC(1) = GUI.APP(1,6,uidropdown(obj.app.maingrid,'BackgroundColor','w','Items',{},'ValueChangedFcn',@obj.cb_slider));
            obj.app.dropC(2) = GUI.APP(1,6,uidropdown(obj.app.maingrid,'BackgroundColor','w','Items',{},'ValueChangedFcn',@obj.cb_slider,'Visible','off'));
            obj.app.grid(3)  = GUI.APP(2,[5 7],uigridlayout(obj.app.maingrid,'RowHeight',{'1x',40,30},'ColumnWidth',{20,150,'1x','1x',120,30,20},'Padding',[15 10 15 0],'RowSpacing',5,'BackgroundColor','w'));
            obj.app.axes     = GUI.APP(1,[2 6],uiaxes(obj.app.grid(3),'BackgroundColor','w','Box','on'));
            obj.app.waittip  = GUI.APP(1,[2 6],uilabel(obj.app.grid(3),'HorizontalAlignment','center','Text','                 Please wait ... ...','Visible','off'));
            tempTb = axtoolbar(obj.app.axes(1),{'rotate','pan','zoomin','zoomout'});
            obj.app.toolC(1)   = axtoolbarbtn(tempTb,'push','Icon',obj.GUI.icon.gif,'Tooltip','Save the evolutionary process to gif','ButtonPushedFcn',@obj.cb_toolbutton1);
            obj.app.toolC(2)   = axtoolbarbtn(tempTb,'push','Icon',obj.GUI.icon.newfigure,'Tooltip','Open in new figure and save to workspace','ButtonPushedFcn',@obj.cb_toolbutton2);
            obj.app.slider     = GUI.APP(2,[1 7],uislider(obj.app.grid(3),'Limits',[0 1],'MajorTicks',0:0.25:1,'MajorTickLabels',{'0%','25%','50%','75%','100%'},'MinorTicks',0:0.01:1,'ValueChangedFcn',@obj.cb_slider));
            obj.app.labelC     = GUI.APP(3,[1 2],uilabel(obj.app.grid(3),'Text','','HorizontalAlignment','left'));
            obj.app.buttonC(1) = GUI.APP(3,3,uibutton(obj.app.grid(3),'push','Text','Start','FontSize',16,'ButtonpushedFcn',@obj.cb_start));
            obj.app.buttonC(2) = GUI.APP(3,4,uibutton(obj.app.grid(3),'push','Text','Stop','FontSize',16,'Enable','off','ButtonpushedFcn',@obj.cb_stop));
            obj.app.menuC      = uicontext(obj.GUI.app.figure,120);
            obj.app.menuC.add('  Save best solutions','',{@obj.cb_save,1});
            obj.app.menuC.add('  Save all solutions','',{@obj.cb_save,2});
            obj.app.menuC.flush();
            obj.app.buttonC(3) = GUI.APP(3,[6 7],uibutton(obj.app.grid(3),'push','Text','Save','FontSize',16,'Enable','off','ButtonpushedFcn',@(~,~)obj.app.menuC.show()));
            
            % The fourth panel
            obj.app.grid(4)  = GUI.APP(2,9,uigridlayout(obj.app.maingrid,'RowHeight',{22,22,'1x'},'ColumnWidth',{'1x','1x'},'Padding',[12 10 12 0],'RowSpacing',15,'BackgroundColor','w'));
            obj.app.dropD(1) = GUI.APP(1,[1 2],uidropdown(obj.app.grid(4),'BackgroundColor','w','Items',{},'ItemsData',1:999,'ValueChangedFcn',@obj.cb_dropdown1));
            obj.app.dropD(2) = GUI.APP(2,1,uidropdown(obj.app.grid(4),'BackgroundColor','w','Items',{},'ValueChangedFcn',@obj.cb_dropdown2));
            obj.app.dropD(3) = GUI.APP(2,1,uidropdown(obj.app.grid(4),'BackgroundColor','w','Items',{},'ValueChangedFcn',@obj.cb_dropdown2,'Visible','off'));
            obj.app.labelD   = GUI.APP(2,2,uilabel(obj.app.grid(4),'Text','','HorizontalAlignment','right'));
            obj.app.textD    = GUI.APP(3,[1 2],uitextarea(obj.app.grid(4),'Editable','off'));
            
            % Initialization
            obj.cb_filter([],[],2);
            obj.app.listA(1).Value = 'NSGAII';
            obj.app.listA(2).Value = 'DTLZ2';
            obj.cb_updateList([],[],1);
            obj.cb_updateList([],[],-1);
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
            show = cellfun(func,obj.GUI.algList(:,1));
            obj.app.listA(1).Items = ['(Open File)';obj.GUI.algList(show,2)];
            obj.app.listA(1).Value = {};
            obj.app.labelA(5).Text = sprintf('%d / %d',sum(show),length(show));
            % Update the list of problems
            show = cellfun(func,obj.GUI.proList(:,1));
            obj.app.listA(2).Items = ['(Open File)';obj.GUI.proList(show,2)];
            obj.app.listA(2).Value = {};
            obj.app.labelA(7).Text = sprintf('%d / %d',sum(show),length(show));
            % Update the lists of metrics
            show   = cellfun(@(s)func(s(2:end,1:end-2)),obj.GUI.metList(:,1));
            if obj.app.stateA(1).Value == 0 % Multi-objective optimization
                obj.app.dropC(1).Items = ['Population (objectives)';'Population (variables)';'True Pareto front';obj.GUI.metList(show,2)];
                obj.app.dropD(2).Items = ['runtime';obj.GUI.metList(show,2)];
                obj.cb_dropdown2();
            else                            % Single-objective optimization
                obj.app.dropC(2).Items = ['Population (variables)';obj.GUI.metList(show,2)];
                obj.app.dropD(3).Items = ['runtime';obj.GUI.metList(show,2)];
                obj.cb_dropdown2();
            end
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
            obj.app.listB.del([],type);
            obj.app.listB.add(filename,type);
            obj.app.listB.flush();
        end
        %% Start the execution
        function cb_start(obj,~,~)
            if strcmp(obj.app.buttonC(1).Text,'Pause')
                obj.app.buttonC(1).Text = 'Continue';
            elseif strcmp(obj.app.buttonC(1).Text,'Continue')
                obj.app.buttonC(1).Text = 'Pause';
            else
                try
                    % Generate the ALGORITHM and PROBLEM object
                    for i = 1 : 2
                        item = obj.app.listB.items(i);
                        para = cell(1,length(item.edit));
                        for j = 1 : length(para)
                            if ~isempty(item.edit(j).Value)
                                para{j} = str2num(item.edit(j).Value);
                                assert(~isempty(para{j}),'The parameter "%s" of %s is illegal.',item.label(j).Text,item.title.Text);
                            end
                        end
                        if i == 1
                            ALG = feval(item.title.Text,'parameter',para,'outputFcn',@obj.outputFcn,'save',inf);
                        else
                            PRO = feval(item.title.Text,'N',para{1},'M',para{2},'D',para{3},item.label(4).Text,para{4},'parameter',para(5:end));
                        end
                    end
                catch err
                    uialert(obj.GUI.app.figure,err.message,'Invalid parameter settings');
                    return;
                end
                % Update the data
                str = sprintf('<Algorithm: %s>\n',class(ALG));
                for i = 1 : length(obj.app.listB.items(1).label)
                    str = [str,sprintf('%s: %s\n',obj.app.listB.items(1).label(i).Text,obj.app.listB.items(1).edit(i).Value)];
                end
                str = [str,sprintf('\n<Problem: %s>\n',class(PRO))];
                for i = 1 : length(obj.app.listB.items(2).label)
                    str = [str,sprintf('%s: %s\n',obj.app.listB.items(2).label(i).Text,obj.app.listB.items(2).edit(i).Value)];
                end
                obj.data = [obj.data;{ALG},{PRO},{str}];
                % Update the GUI
                [obj.GUI.app.button.Enable]  = deal('off');
                [obj.app.stateA.Enable]      = deal('off');
                [obj.app.listA.Enable]       = deal('off');
                [obj.app.toolC(1:2).Visible] = deal('off');
                [obj.app.dropD.Enable]       = deal('off');
                obj.app.listB.Enable         = 'off';
                obj.app.buttonC(1).Text      = 'Pause';
                obj.app.buttonC(2).Enable    = 'on';
                obj.app.buttonC(3).Enable    = 'on';
                if PRO.M > 1
                    obj.app.dropC(1).Value   = obj.app.dropC(1).Items{1};
                    obj.app.dropC(1).Visible = 'on';
                    obj.app.dropC(2).Visible = 'off';
                    obj.app.dropD(2).Visible = 'on';
                    obj.app.dropD(3).Visible = 'off';
                else
                    obj.app.dropC(2).Value   = obj.app.dropC(2).Items{1};
                    obj.app.dropC(1).Visible = 'off';
                    obj.app.dropC(2).Visible = 'on';
                    obj.app.dropD(2).Visible = 'off';
                    obj.app.dropD(3).Visible = 'on';
                end
                obj.app.dropC(1).Enable = 'off';
                obj.app.dropC(2).Enable = 'off';
                obj.app.dropD(1).Items  = [obj.app.dropD(1).Items,sprintf('%s on %s',class(ALG),class(PRO))];
                obj.app.dropD(1).Value  = length(obj.app.dropD(1).Items);
                obj.app.labelD.Text     = '';
                obj.app.textD.Value     = '';
                % Execute the algorithm
                try
                    ALG.Solve(PRO);
                    obj.cb_stop();
                catch err
                    uialert(obj.GUI.app.figure,'The algorithm terminates unexpectedly, please refer to the command window for details.','Error');
                    obj.cb_stop();
                    rethrow(err);
                end
            end
        end
        %% Stop the execution
        function cb_stop(obj,~,~)
            [obj.GUI.app.button.Enable]  = deal('on');
            [obj.app.stateA.Enable]      = deal('on');
            [obj.app.listA.Enable]       = deal('on');
            [obj.app.toolC(1:2).Visible] = deal('on');
            [obj.app.dropD.Enable]       = deal('on');
            obj.app.listB.Enable         = 'on';
            obj.app.buttonC(1).Text      = 'Start';
            obj.app.buttonC(2).Enable    = 'off';
            obj.app.dropC(1).Enable      = 'on';
            obj.app.dropC(2).Enable      = 'on';
            if isempty(obj.data{end,1}.result)
                obj.data(end,:)             = [];
                obj.app.dropD(1).Items(end) = [];
            end
            obj.cb_dropdown1();
        end
        %% Save the result
        function cb_save(obj,~,~,type)
            ALG   = obj.data{obj.app.dropD(1).Value,1};
            PRO   = obj.data{obj.app.dropD(1).Value,2};
            rate  = PRO.FE/max(PRO.FE,PRO.maxFE);
            index = max(1,round(obj.app.slider.Value/rate*size(ALG.result,1)));
            Pop   = ALG.result{index,2};
            if type == 1
                Pop = Pop.best;
            end
            if isempty(Pop)
                uialert(obj.GUI.app.figure,'No solution can be saved.','Error');
            else
                Data = [Pop.decs,Pop.objs,Pop.cons];
                try
                    [Name,Path] = uiputfile({'*.txt','Text file';'*.dat','Text file';'*.csv','Text file';'*.mat','MAT file';'*.xlsx','Excel table'},'','data');
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
                    uialert(obj.GUI.app.figure,'Fail to save the result, please refer to the command window for details.','Error');
                    rethrow(err);
                end
            end
        end
        %% Output function
        function outputFcn(obj,Algorithm,Problem)
            obj.app.slider.Value = Problem.FE/max(Problem.FE,Problem.maxFE);
            obj.cb_slider();
            assert(strcmp(obj.app.buttonC(2).Enable,'on'),'PlatEMO:Termination','');
            if strcmp(obj.app.buttonC(1).Text,'Continue')
                waitfor(obj.app.buttonC(1),'Text');
            end
            assert(strcmp(obj.app.buttonC(2).Enable,'on'),'PlatEMO:Termination','');
        end
        %% Show the specified data
        function cb_slider(obj,~,~,ax)
            if ~isempty(obj.app.dropD(1).Items)
                % Determine the current number of evaluationsnumber of evaluations
                ALG  = obj.data{obj.app.dropD(1).Value,1};
                PRO  = obj.data{obj.app.dropD(1).Value,2};
                rate = PRO.FE/max(PRO.FE,PRO.maxFE);
                obj.app.slider.Value      = min(obj.app.slider.Value,rate);
                obj.app.slider.MajorTicks = 0 : 0.25 : rate;
                obj.app.slider.MinorTicks = 0 : 0.01 : rate;
                index = max(1,round(obj.app.slider.Value/rate*size(ALG.result,1)));
                obj.app.labelC.Text = sprintf('%d evaluations',ALG.result{index,1});
                % Clear the default or specified axes
                if nargin > 3
                    Draw(ax);
                else
                    Draw(obj.app.axes);
                end
                isMetric = false;
                if PRO.M > 1    % Multi-objective optimization
                    switch obj.app.dropC(1).Value
                        case 'Population (objectives)'
                            PRO.DrawObj(ALG.result{index,2});
                        case 'Population (variables)'
                            PRO.DrawDec(ALG.result{index,2});
                        case 'True Pareto front'
                            Draw(PRO.optimum,{'\it f\rm_1','\it f\rm_2','\it f\rm_3'});
                        otherwise
                            obj.app.waittip.Visible = 'on'; drawnow();
                            Draw([cell2mat(ALG.result(:,1)),ALG.CalMetric(obj.app.dropC(1).Value)],'-k.','LineWidth',1.5,'MarkerSize',10,{'Number of function evaluations',strrep(obj.app.dropC(1).Value,'_',' '),[]});
                            obj.app.waittip.Visible = 'off';
                            isMetric = true;
                    end
                    if ~isMetric
                        obj.app.dropD(2).Value = 'runtime';
                        obj.app.labelD.Text = sprintf('%.4fs',ALG.CalMetric('runtime'));
                    else
                        obj.app.dropD(2).Value = obj.app.dropC(1).Value;
                        value = ALG.CalMetric(obj.app.dropD(2).Value);
                        obj.app.labelD.Text = sprintf('%.4e',value(end));
                    end
                else            % Single-objective optimization
                    switch obj.app.dropC(2).Value
                        case 'Population (variables)'
                            PRO.DrawDec(ALG.result{index,2});
                        otherwise
                            Draw([cell2mat(ALG.result(:,1)),ALG.CalMetric(obj.app.dropC(2).Value)],'-k.','LineWidth',1.5,'MarkerSize',10,{'Number of function evaluations',strrep(obj.app.dropC(2).Value,'_',' '),[]});
                            isMetric = true;
                    end
                    if ~isMetric
                        obj.app.dropD(3).Value = 'runtime';
                        obj.app.labelD.Text = sprintf('%.4fs',ALG.CalMetric('runtime'));
                    else
                        obj.app.dropD(3).Value = obj.app.dropC(2).Value;
                        value = ALG.CalMetric(obj.app.dropD(3).Value);
                        obj.app.labelD.Text = sprintf('%.4e',value(end));
                    end
                end
            end
        end
        %% Create the gif
        function cb_toolbutton1(obj,~,~)
            if ~isempty(obj.app.dropD(1).Items)
                [file,folder] = uiputfile({'*.gif','GIF image'},'');
                if file ~= 0
                    try
                        filename = fullfile(folder,file);
                        figure('NumberTitle','off','Name','Figure for creating the gif');
                        for i = linspace(0,1,20)
                            obj.app.slider.Value = i;
                            obj.cb_slider([],[],gca);
                            drawnow('limitrate');
                            [I,map] = rgb2ind(print('-RGBImage'),20);
                            if i == 0
                                imwrite(I,map,filename,'gif','Loopcount',inf,'DelayTime',0.2);
                            else
                                imwrite(I,map,filename,'gif','WriteMode','append','DelayTime',0.2);
                            end
                        end
                        delete(gcf);
                    catch
                        uialert(obj.GUI.app.figure,sprintf('Fail to save the gif to %s.',filename),'Error');
                        return;
                    end
                end
            end
        end
        %% Open in new figure
        function cb_toolbutton2(obj,~,~)
            if ~isempty(obj.app.dropD(1).Items)
                if isempty(get(gcf,'CurrentAxes'))
                    axes('FontName',obj.app.axes.FontName,'FontSize',obj.app.axes.FontSize,'NextPlot',obj.app.axes.NextPlot,'Box',obj.app.axes.Box,'View',obj.app.axes.View);
                    copyobj(obj.app.axes.Children,gca);
                else
                    h = copyobj(obj.app.axes.Children,gca);
                    for i = 1 : length(h)
                        if strcmp(h(i).Type,'line')
                        	set(h(i),'Color',rand(1,3),'Markerfacecolor',rand(1,3));
                        end
                    end
                end
                axis tight;
                Data = arrayfun(@(s){s.XData,s.YData,s.ZData},get(gca,'Children'),'UniformOutput',false);
                assignin('base','Data',cat(1,Data{:}));
            end
        end
        %% Change the selected result
        function cb_dropdown1(obj,~,~)
            if ~isempty(obj.app.dropD(1).Items)
                obj.cb_slider();
                obj.app.textD.Value = obj.data{obj.app.dropD(1).Value,3};
                if obj.data{obj.app.dropD(1).Value,2}.M > 1
                    obj.app.dropC(1).Visible = 'on';
                    obj.app.dropC(2).Visible = 'off';
                    obj.app.dropD(2).Visible = 'on';
                    obj.app.dropD(3).Visible = 'off';
                else
                    obj.app.dropC(1).Visible = 'off';
                    obj.app.dropC(2).Visible = 'on';
                    obj.app.dropD(2).Visible = 'off';
                    obj.app.dropD(3).Visible = 'on';
                end
            end
        end
        %% Change the selected metric
        function cb_dropdown2(obj,~,~)
            if ~isempty(obj.app.dropD(1).Items)
                if obj.data{obj.app.dropD(1).Value,2}.M > 1
                    if strcmp(obj.app.dropD(2).Value,'runtime')
                        obj.app.dropC(1).Value = obj.app.dropC(1).Items{1};
                    elseif ~strcmp(obj.app.dropD(2).Value,'runtime')
                        obj.app.dropC(1).Value = obj.app.dropD(2).Value;
                    end
                else
                    if strcmp(obj.app.dropD(3).Value,'runtime')
                        obj.app.dropC(2).Value = obj.app.dropC(2).Items{1};
                    elseif ~strcmp(obj.app.dropD(3).Value,'runtime')
                        obj.app.dropC(2).Value = obj.app.dropD(3).Value;
                    end
                end
                obj.cb_slider();
            end
        end
    end
end