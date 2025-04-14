classdef module_test < handle
%module_test - Test module.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
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
            obj.app.grid(1)   = GUI.APP(2,1,uigridlayout(obj.app.maingrid,'RowHeight',{16,21,16,21,21,16,21,21,21,4,18,'1x',4,18,'1x'},'ColumnWidth',{'1x','1.1x','1x'},'Padding',[8 10 8 0],'RowSpacing',3,'ColumnSpacing',5,'BackgroundColor','w'));
            [obj.app.stateA,obj.app.labelA] = GUI.GenerateLabelButton(obj.app.grid(1),[0,1,0,1,zeros(1,13)],@obj.cb_filter);
            obj.app.labelA(4) = GUI.APP(11,[1 2],uilabel(obj.app.grid(1),'Text','Algorithms','FontSize',13,'FontColor',[.2 .4 .7],'FontWeight','bold'));
            obj.app.labelA(5) = GUI.APP(11,3,uilabel(obj.app.grid(1),'HorizontalAlignment','right','FontSize',10,'FontColor',[.2 .4 .7]));
            obj.app.listA(1)  = GUI.APP(12,[1 3],uilistbox(obj.app.grid(1),'FontColor',[.2 .4 .7]));
            obj.app.labelA(6) = GUI.APP(14,[1 2],uilabel(obj.app.grid(1),'Text','Problems','FontSize',13,'FontColor',[.9 .5 .2],'FontWeight','bold'));
            obj.app.labelA(7) = GUI.APP(14,3,uilabel(obj.app.grid(1),'HorizontalAlignment','right','FontSize',10,'FontColor',[.9 .5 .2]));
            obj.app.listA(2)  = GUI.APP(15,[1 3],uilistbox(obj.app.grid(1),'FontColor',[.9 .5 .2]));
            obj.app.dropA(1)  = GUI.APP(11,2,uidropdown(obj.app.grid(1),'BackgroundColor','w','FontColor',[.2 .4 .7],'Items',{'All year'},'ValueChangedFcn',@(h,~)GUI.UpdateAlgProListYear(obj.app.listA(1),h,obj.app.labelA(5),obj.GUI.algList)));
            obj.app.dropA(2)  = GUI.APP(14,2,uidropdown(obj.app.grid(1),'BackgroundColor','w','FontColor',[.9 .5 .2],'Items',{'All year'},'ValueChangedFcn',@(h,~)GUI.UpdateAlgProListYear(obj.app.listA(2),h,obj.app.labelA(7),obj.GUI.proList)));

            % The second panel
            obj.app.listB   = uilist(obj.app.maingrid,obj.GUI.app.figure,obj.GUI.icon);
            obj.app.grid(2) = GUI.APP(2,3,obj.app.listB.grid);
            obj.app.listA(1).ValueChangedFcn = @(~,~)GUI.UpdateAlgProPara(obj.GUI.app.figure,obj.app.listA(1),obj.app.listB,'ALGORITHM',1);
            obj.app.listA(2).ValueChangedFcn = @(~,~)GUI.UpdateAlgProPara(obj.GUI.app.figure,obj.app.listA(2),obj.app.listB,'PROBLEM',-1);

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
            obj.app.buttonC(2) = GUI.APP(3,4,uibutton(obj.app.grid(3),'push','Text','Stop','FontSize',16,'Enable','off','ButtonpushedFcn',@(~,~)set(obj.app.buttonC(1:2),{'Enable','Text'},{true,'Start';false,'Stop'})));
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
            GUI.UpdateAlgProPara(obj.GUI.app.figure,obj.app.listA(1),obj.app.listB,'ALGORITHM',1);
            GUI.UpdateAlgProPara(obj.GUI.app.figure,obj.app.listA(2),obj.app.listB,'PROBLEM',-1);
        end
    end
    methods(Access = private)
        %% Update the algorithms and problems in the lists
        function cb_filter(obj,~,~,index)
            % Update the lists of algorithms and problems
            func = GUI.UpdateAlgProList(index,obj.app.stateA,obj.app.listA(1),obj.app.dropA(1),obj.app.labelA(5),obj.GUI.algList,obj.app.listA(2),obj.app.dropA(2),obj.app.labelA(7),obj.GUI.proList);
            % Update the list of metrics
            show = cellfun(@(s)func(s(2:end,1:end-2)),obj.GUI.metList(:,1));
            if obj.app.stateA(1).Value == 0 % Multi-objective optimization
                obj.app.dropC(1).Items = ['Population (objectives)';'Population (variables)';'True Pareto front';obj.GUI.metList(show,2)];
                obj.app.dropD(2).Items = ['runtime';obj.GUI.metList(show,2)];
            else                            % Single-objective optimization
                obj.app.dropC(2).Items = ['Population (variables)';obj.GUI.metList(show,2)];
                obj.app.dropD(3).Items = ['runtime';obj.GUI.metList(show,2)];
            end
            obj.cb_dropdown2();
        end
        %% Start the execution
        function cb_start(obj,~,~)
            if strcmp(obj.app.buttonC(1).Text,'Pause')
                obj.app.buttonC(1).Text = 'Continue';
            elseif strcmp(obj.app.buttonC(1).Text,'Continue')
                obj.app.buttonC(1).Text = 'Pause';
            else
                % Generate the ALGORITHM and PROBLEM objects
                try
                    [name,para] = GUI.GetParameterSetting(obj.app.listB.items(1));
                    ALG = feval(name,'parameter',para,'outputFcn',@obj.outputFcn,'save',20);
                    [name,para] = GUI.GetParameterSetting(obj.app.listB.items(2));
                    PRO = feval(name,'N',para{1},'M',para{2},'D',para{3},obj.app.listB.items(2).label(4).Text,para{4},'parameter',para(5:end));
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
                set([obj.GUI.app.button,obj.app.stateA,obj.app.listA,obj.app.dropA,obj.app.dropD,obj.app.dropC],'Enable',false);
                obj.app.listB.Enable = false;
                set(obj.app.toolC,'Visible',false);
                obj.app.buttonC(1).Text = 'Pause';
                set(obj.app.buttonC([2,3]),'Enable',true);
                if PRO.M > 1
                    obj.app.dropC(1).Value = obj.app.dropC(1).Items{1};
                    set([obj.app.dropC(1),obj.app.dropD(2)],'Visible',true);
                    set([obj.app.dropC(2),obj.app.dropD(3)],'Visible',false);
                else
                    obj.app.dropC(2).Value = obj.app.dropC(2).Items{1};
                    set([obj.app.dropC(1),obj.app.dropD(2)],'Visible',false);
                    set([obj.app.dropC(2),obj.app.dropD(3)],'Visible',true);
                end
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
            set([obj.GUI.app.button,obj.app.stateA,obj.app.listA,obj.app.dropA,obj.app.dropD,obj.app.dropC],'Enable',true);
            obj.app.listB.Enable = true;
            set(obj.app.toolC,'Visible',true);
            obj.app.buttonC(1).Text   = 'Start';
            obj.app.buttonC(2).Enable = false;
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
            GUI.SavePopulation(obj.GUI.app.figure,ALG.result{index,2},type);
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
                % Determine the current number of evaluations
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
                figure(obj.GUI.app.figure);
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
                    set([obj.app.dropC(1),obj.app.dropD(2)],'Visible',true);
                    set([obj.app.dropC(2),obj.app.dropD(3)],'Visible',false);
                else
                    set([obj.app.dropC(1),obj.app.dropD(2)],'Visible',false);
                    set([obj.app.dropC(2),obj.app.dropD(3)],'Visible',true);
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