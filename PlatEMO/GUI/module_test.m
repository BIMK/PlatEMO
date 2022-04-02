classdef module_test < handle
%module_test - Test module.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
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
            obj.app.maingrid = GUI.APP(3,1,uigridlayout(obj.GUI.app.maingrid,'RowHeight',{20,'1x'},'ColumnWidth',{'1.2x',1,'1x',1,'3x',1,'1x'},'Padding',[0 5 0 5],'RowSpacing',5,'ColumnSpacing',0,'BackgroundColor','w'));
            obj.app.label(1) = GUI.APP(1,1,uilabel(obj.app.maingrid,'Text','Algorithm selection','HorizontalAlignment','center','FontSize',11,'FontColor',[.5 .5 .5]));
            obj.app.label(2) = GUI.APP(1,3,uilabel(obj.app.maingrid,'Text','Parameter setting','HorizontalAlignment','center','FontSize',11,'FontColor',[.5 .5 .5]));
            obj.app.label(3) = GUI.APP(1,5,uilabel(obj.app.maingrid,'Text','Result display','HorizontalAlignment','center','FontSize',11,'FontColor',[.5 .5 .5]));
            obj.app.label(4) = GUI.APP(1,7,uilabel(obj.app.maingrid,'Text','Result selection','HorizontalAlignment','center','FontSize',11,'FontColor',[.5 .5 .5]));
            GUI.APP([1 2],2,uipanel(obj.app.maingrid,'BackgroundColor',[.8 .8 .8]));
            GUI.APP([1 2],4,uipanel(obj.app.maingrid,'BackgroundColor',[.8 .8 .8]));
            GUI.APP([1 2],6,uipanel(obj.app.maingrid,'BackgroundColor',[.8 .8 .8]));

            % The first panel
            obj.app.grid(1)    = GUI.APP(2,1,uigridlayout(obj.app.maingrid,'RowHeight',{16,20,16,20,16,20,20,25,'1x',25,'1x'},'ColumnWidth',{'1x','1x','1x'},'Padding',[8 10 8 0],'RowSpacing',5,'ColumnSpacing',5,'BackgroundColor','w'));
            obj.app.labelA(1)  = GUI.APP(1,[1 3],uilabel(obj.app.grid(1),'Text','Number of objectives','VerticalAlignment','bottom','FontSize',12,'FontColor',[.15 .6 .2],'FontWeight','bold'));
            obj.app.stateA(1)  = GUI.APP(2,1,uibutton(obj.app.grid(1),'state','Text','single','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',0,'Tooltip','The problem has a single objective','ValueChangedFcn',{@obj.cb_filter,1}));
            obj.app.stateA(2)  = GUI.APP(2,2,uibutton(obj.app.grid(1),'state','Text','multi','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',1,'Tooltip','The problem has 2 or 3 objectives','ValueChangedFcn',{@obj.cb_filter,2}));
            obj.app.stateA(3)  = GUI.APP(2,3,uibutton(obj.app.grid(1),'state','Text','many','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',0,'Tooltip','The problem has more than 3 objectives','ValueChangedFcn',{@obj.cb_filter,3}));
            obj.app.labelA(2)  = GUI.APP(3,[1 3],uilabel(obj.app.grid(1),'Text','Encoding scheme','VerticalAlignment','bottom','FontSize',12,'FontColor',[.15 .6 .2],'FontWeight','bold'));
            obj.app.stateA(4)  = GUI.APP(4,1,uibutton(obj.app.grid(1),'state','Text','real','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',1,'Tooltip','The decision variables are real values','ValueChangedFcn',{@obj.cb_filter,4}));
            obj.app.stateA(5)  = GUI.APP(4,2,uibutton(obj.app.grid(1),'state','Text','binary','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',0,'Tooltip','The decision variables are binary values','ValueChangedFcn',{@obj.cb_filter,5}));
            obj.app.stateA(6)  = GUI.APP(4,3,uibutton(obj.app.grid(1),'state','Text','permutation','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',0,'Tooltip','The decision vector is a permutation','ValueChangedFcn',{@obj.cb_filter,6}));
            obj.app.labelA(3)  = GUI.APP(5,[1 3],uilabel(obj.app.grid(1),'Text','Special difficulties','VerticalAlignment','bottom','FontSize',12,'FontColor',[.15 .6 .2],'FontWeight','bold'));
            obj.app.stateA(7)  = GUI.APP(6,1,uibutton(obj.app.grid(1),'state','Text','large','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',0,'Tooltip','The problem has more than 100 decision variables','ValueChangedFcn',{@obj.cb_filter,7}));
            obj.app.stateA(8)  = GUI.APP(6,2,uibutton(obj.app.grid(1),'state','Text','constrained','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',0,'Tooltip','The problem has constraints','ValueChangedFcn',{@obj.cb_filter,8}));
            obj.app.stateA(9)  = GUI.APP(6,3,uibutton(obj.app.grid(1),'state','Text','expensive','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',0,'Tooltip','The objectives are computationally time-consuming','ValueChangedFcn',{@obj.cb_filter,9}));
            obj.app.stateA(10) = GUI.APP(7,1,uibutton(obj.app.grid(1),'state','Text','multimodal','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',0,'Tooltip','The objectives are multimodal','ValueChangedFcn',{@obj.cb_filter,10}));
            obj.app.stateA(11) = GUI.APP(7,2,uibutton(obj.app.grid(1),'state','Text','sparse','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',0,'Tooltip','Most decision variables of the optimal solutions are zero','ValueChangedFcn',{@obj.cb_filter,11}));
            obj.app.stateA(12) = GUI.APP(7,3,uibutton(obj.app.grid(1),'state','Text','preference','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Value',0,'Tooltip','Only searching for preferred regions of the Pareto front','ValueChangedFcn',{@obj.cb_filter,12}));
            obj.app.labelA(4)  = GUI.APP(8,[1 2],uilabel(obj.app.grid(1),'Text','Algorithms','VerticalAlignment','bottom','FontSize',15,'FontColor',[.2 .4 .7],'FontWeight','bold'));
            obj.app.labelA(5)  = GUI.APP(8,3,uilabel(obj.app.grid(1),'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',10,'FontColor',[.2 .4 .7]));
            obj.app.listA(1)   = GUI.APP(9,[1 3],uilistbox(obj.app.grid(1),'FontColor',[.2 .4 .7],'ValueChangedFcn',{@obj.cb_updateList,1}));
            obj.app.labelA(6)  = GUI.APP(10,[1 2],uilabel(obj.app.grid(1),'Text','Problems','VerticalAlignment','bottom','FontSize',15,'FontColor',[.9 .5 .2],'FontWeight','bold'));
            obj.app.labelA(7)  = GUI.APP(10,3,uilabel(obj.app.grid(1),'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',10,'FontColor',[.9 .5 .2]));
            obj.app.listA(2)   = GUI.APP(11,[1 3],uilistbox(obj.app.grid(1),'FontColor',[.9 .5 .2],'ValueChangedFcn',{@obj.cb_updateList,-1}));

            % The second panel
            obj.app.listB   = uilist(obj.app.maingrid,obj.GUI.app.figure,obj.GUI.icon);
            obj.app.grid(2) = GUI.APP(2,3,obj.app.listB.grid);

            % The third panel
            obj.app.grid(3) = GUI.APP(2,5,uigridlayout(obj.app.maingrid,'RowHeight',{'1x',40,30},'ColumnWidth',{20,150,'1x','1x',150,20},'Padding',[15 10 15 0],'RowSpacing',5,'BackgroundColor','w'));
            obj.app.axes    = GUI.APP(1,[2 5],uiaxes(obj.app.grid(3),'BackgroundColor','w','Box','on'));
            obj.app.waittip = GUI.APP(1,[2 5],uilabel(obj.app.grid(3),'HorizontalAlignment','center','Text','                 Please wait ... ...','Visible','off'));
            obj.app.menu(1) = uicontext2(obj.GUI.app.figure,@obj.cb_slider);
            obj.app.menu(1).add('Population (obj.)',false);
            obj.app.menu(1).add('Population (dec.)',false);
            obj.app.menu(1).add('True Pareto front',true);
            for i = 1 : size(obj.GUI.metList,1)
                obj.app.menu(1).add(obj.GUI.metList{i,2},false);
            end
            obj.app.menu(1).add('Feasible rate',false);
            obj.app.menu(1).flush();
            obj.app.menu(2) = uicontext2(obj.GUI.app.figure,@obj.cb_slider);
            obj.app.menu(2).add('Population (dec.)',true);
            obj.app.menu(2).add('Minimum value',false);
            obj.app.menu(2).add('Feasible rate',false);
            obj.app.menu(2).flush();
            tempTb = axtoolbar(obj.app.axes(1),{'rotate','pan','zoomin','zoomout'});
            obj.app.toolC(1)   = axtoolbarbtn(tempTb,'push','Icon',obj.GUI.icon.gif,'Tooltip','Save the evolutionary process to gif','ButtonPushedFcn',@obj.cb_toolbutton1);
            obj.app.toolC(2)   = axtoolbarbtn(tempTb,'push','Icon',obj.GUI.icon.newfigure,'Tooltip','Open in new figure and save to workspace','ButtonPushedFcn',@obj.cb_toolbutton2);
            obj.app.toolC(3)   = axtoolbarbtn(tempTb,'push','Icon',obj.GUI.icon.datasource,'Tooltip','Data source','ButtonPushedFcn',@obj.cb_toolbutton3);
            obj.app.slider     = GUI.APP(2,[1 6],uislider(obj.app.grid(3),'Limits',[0 1],'MajorTicks',0:0.25:1,'MajorTickLabels',{'0%','25%','50%','75%','100%'},'MinorTicks',0:0.01:1,'ValueChangedFcn',@obj.cb_slider));
            obj.app.labelC     = GUI.APP(3,[5 6],uilabel(obj.app.grid(3),'Text','','HorizontalAlignment','right'));
            obj.app.buttonC(1) = GUI.APP(3,3,uibutton(obj.app.grid(3),'push','Text','Start','FontSize',16,'ButtonpushedFcn',@obj.cb_start));
            obj.app.buttonC(2) = GUI.APP(3,4,uibutton(obj.app.grid(3),'push','Text','Stop','FontSize',16,'Enable','off','ButtonpushedFcn',@obj.cb_stop));
            
            % The fourth panel
            obj.app.grid(4)  = GUI.APP(2,7,uigridlayout(obj.app.maingrid,'RowHeight',{22,22,'1x'},'ColumnWidth',{'1x','1x'},'Padding',[12 10 12 0],'RowSpacing',15,'BackgroundColor','w'));
            obj.app.dropD(1) = GUI.APP(1,[1 2],uidropdown(obj.app.grid(4),'BackgroundColor','w','Items',{},'ItemsData',1:999,'ValueChangedFcn',@obj.cb_dropdown1));
            obj.app.dropD(2) = GUI.APP(2,1,uidropdown(obj.app.grid(4),'BackgroundColor','w','Items',['runtime';obj.GUI.metList(:,2);'Feasible rate'],'ValueChangedFcn',@obj.cb_dropdown2));
            obj.app.dropD(3) = GUI.APP(2,1,uidropdown(obj.app.grid(4),'BackgroundColor','w','Items',{'runtime';'Minimum value';'Feasible rate'},'ValueChangedFcn',@obj.cb_dropdown2,'Visible','off'));
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
            elseif index < 7
                [obj.app.stateA(4:6).Value] = deal(0);
                obj.app.stateA(index).Value = 1;
            end
            filter = [obj.app.stateA.Value];
            func   = @(s)all(any(repmat([true,filter],size(s,1),1)&s,2)) && all((any(s(:,2:end),1)&filter)==filter);
            show   = cellfun(func,obj.GUI.algList(:,1));
            obj.app.listA(1).Items = ['(Open File)';obj.GUI.algList(show,2)];
            obj.app.listA(1).Value = {};
            obj.app.labelA(5).Text = sprintf('%d / %d',sum(show),length(show));
            show   = cellfun(func,obj.GUI.proList(:,1));
            obj.app.listA(2).Items = ['(Open File)';obj.GUI.proList(show,2)];
            obj.app.listA(2).Value = {};
            obj.app.labelA(7).Text = sprintf('%d / %d',sum(show),length(show));
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
                [file,path] = uigetfile('*.m','');
                if file ~= 0
                    try
                        filename = fullfile(path,file);
                        f   = fopen(filename);
                        str = fgetl(f);
                        fclose(f);
                        assert(contains(str,['< ',tip]));
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
                            PRO = feval(item.title.Text,'N',para{1},'M',para{2},'D',para{3},'maxFE',para{4},'parameter',para(5:end));
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
                obj.app.menu((PRO.M<=1)+1).value = 1;
                [obj.app.menu(1).items(3:end).Enable] = deal('off');
                [obj.app.menu(2).items(2:end).Enable] = deal('off');
                obj.app.dropD(1).Items = [obj.app.dropD(1).Items,sprintf('%s on %s',class(ALG),class(PRO))];
                obj.app.dropD(1).Value = length(obj.app.dropD(1).Items);
                obj.app.labelD.Text    = '';
                obj.app.textD.Value    = '';
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
            [obj.app.menu(1).items(3:end).Enable] = deal('on');
            [obj.app.menu(2).items(2:end).Enable] = deal('on');
            if isempty(obj.data{end,1}.result)
                obj.data(end,:)             = [];
                obj.app.dropD(1).Items(end) = [];
            end
            obj.cb_dropdown1();
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
                if PRO.M > 1
                    % Show the result of multi-objective optimization
                    switch obj.app.menu(1).value
                        case 1
                            PRO.DrawObj(ALG.result{index,2});
                        case 2
                            PRO.DrawDec(ALG.result{index,2});
                        case 3
                            Draw(PRO.optimum,{'\it f\rm_1','\it f\rm_2','\it f\rm_3'});
                        otherwise
                            obj.app.waittip.Visible = 'on'; drawnow();
                            Draw(ALG.Metric(obj.app.menu(1).string,10),'-k.','LineWidth',1.5,'MarkerSize',10,{'Number of function evaluations',obj.app.menu(1).string,[]});
                            obj.app.waittip.Visible = 'off';
                    end
                    if obj.app.menu(1).value < 4
                        obj.app.dropD(2).Value = 'runtime';
                        obj.app.labelD.Text = sprintf('%.4fs',ALG.Metric('runtime'));
                    else
                        obj.app.dropD(2).Value = obj.app.menu(1).string;
                        value = ALG.Metric(obj.app.dropD(2).Value,10);
                        obj.app.labelD.Text = sprintf('%.4e',value(end));
                    end
                else
                    % Show the result of single-objective optimization
                    switch obj.app.menu(2).value
                        case 1
                            PRO.DrawDec(ALG.result{index,2});
                        otherwise
                            Draw(ALG.Metric(obj.app.menu(2).string,10),'-k.','LineWidth',1.5,'MarkerSize',10,{'Number of function evaluations',obj.app.menu(2).string,[]});
                    end
                    if obj.app.menu(2).value < 2
                        obj.app.dropD(3).Value = 'runtime';
                        obj.app.labelD.Text = sprintf('%.4fs',ALG.Metric('runtime'));
                    else
                        obj.app.dropD(3).Value = obj.app.menu(2).string;
                        value = ALG.Metric(obj.app.dropD(3).Value,10);
                        obj.app.labelD.Text = sprintf('%.4e',value(end));
                    end
                end
            end
        end
        %% Create the gif
        function cb_toolbutton1(obj,~,~)
            if ~isempty(obj.app.dropD(1).Items)
                [file,folder] = uiputfile('*.gif','');
                if file ~= 0
                    try
                        filename = fullfile(folder,file);
                        figure('NumberTitle','off','Name','Figure for creating the gif');
                        for i = linspace(0,1,20)
                            obj.app.slider.Value = i;
                            obj.cb_slider([],[],gca);
                            drawnow();
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
        %% Show the menu of data source
        function cb_toolbutton3(obj,~,~)
            if ~isempty(obj.app.dropD(1).Items)
                obj.app.menu((obj.data{obj.app.dropD(1).Value,2}.M<=1)+1).show();
            end
        end
        %% Change the selected result
        function cb_dropdown1(obj,~,~)
            if ~isempty(obj.app.dropD(1).Items)
                obj.cb_slider();
                obj.app.textD.Value = obj.data{obj.app.dropD(1).Value,3};
                if obj.data{obj.app.dropD(1).Value,2}.M > 1
                    obj.app.dropD(2).Visible = 'on';
                    obj.app.dropD(3).Visible = 'off';
                else
                    obj.app.dropD(2).Visible = 'off';
                    obj.app.dropD(3).Visible = 'on';
                end
            end
        end
        %% Change the selected metric
        function cb_dropdown2(obj,~,~)
            if ~isempty(obj.app.dropD(1).Items)
                if obj.data{obj.app.dropD(1).Value,2}.M > 1
                    if strcmp(obj.app.dropD(2).Value,'runtime') && obj.app.menu(1).value > 3
                        obj.app.menu(1).value = 1;
                    elseif ~strcmp(obj.app.dropD(2).Value,'runtime')
                        obj.app.menu(1).value = find(ismember(obj.app.dropD(2).Items,obj.app.dropD(2).Value),1) + 2;
                    end
                else
                    if strcmp(obj.app.dropD(3).Value,'runtime') && obj.app.menu(2).value > 1
                        obj.app.menu(2).value = 1;
                    elseif ~strcmp(obj.app.dropD(3).Value,'runtime')
                        obj.app.menu(2).value = find(ismember(obj.app.dropD(3).Items,obj.app.dropD(3).Value),1);
                    end
                end
                obj.cb_slider();
            end
        end
    end
end