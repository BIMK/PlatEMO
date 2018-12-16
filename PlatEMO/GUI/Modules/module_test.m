classdef module_test < module
%module_test - Test module.

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
        function obj = module_test(GUI)
            Panel = newMpanel(GUI.figure,[1 1 1200 549,1 1 1 1],4);
            Panel.move(1,-146);
            Panel.move(2,-242);
            Panel.move(3,72);
            obj = obj@module(GUI,Panel);
            
            % The first panel
            panel = Panel.panels(1);
            note = newLabel(panel,[10 10 140 115,0 0 1 1],'','FontSize',9,'HorizontalAlignment','left');
            newLabel(panel,[15 530 135 15,0 0 0 1],'Function selection','ForegroundColor',[.7 .7 .7],'FontSize',9);
            obj.control.setLabel(1) = newLabel(panel,[10 420 135 50,0 0 0 1],'Algorithm','ForegroundColor',[.2 .4 .7],'FontSize',12,'FontWeight','bold','movecallback',@(~,~)set(note.handle,'String','Select an algorithm to be executed'));
            obj.control.setLabel(2) = newLabel(panel,[10 310 135 50,0 0 0 1],'Problem','ForegroundColor',[.9 .5 .2],'FontSize',12,'FontWeight','bold','movecallback',@(~,~)set(note.handle,'String','Select a problem to be solved'));
            obj.control.algPopmenu  = newPopmenu(panel,[20 420 115 20,0 0 0 1],[.2 .4 .7],obj.GUI.algList,'callback',@(~,~)obj.cb_apoPopmenu('alg'));
            obj.control.proPopmenu  = newPopmenu(panel,[20 310 115 20,0 0 0 1],[.9 .5 .2],obj.GUI.proList,'callback',@(~,~)obj.cb_apoPopmenu('pro'));
            obj.control.runButton   = newButton(panel,[40 180 75 40,0 0 0 1],'RUN','FontSize',20,'callback',@obj.cb_run);
            newLine(panel,[5 130 145 1,0 0 0 1]);
            
            % The second panel
            panel = Panel.panels(2);
            newLabel(panel,[15 530 180 15,0 0 0 1],'Parameter setting','ForegroundColor',[.7 .7 .7],'FontSize',9);
            obj.control.setPanel = ParameterList(panel,[6 10 194 514,1 1 1 1],note,obj.GUI.icons);
            obj.control.algPopmenu.index = find(strcmp('NSGAII',get([obj.control.algPopmenu.items.handle],'String')),1);
            obj.control.proPopmenu.index = find(strcmp('DTLZ2',get([obj.control.proPopmenu.items.handle],'String')),1);
            
            % The third panel
            panel = Panel.panels(3);
            newLabel(panel,[255 530 100 15,0 0 0 1],'Result display','ForegroundColor',[.7 .7 .7],'FontSize',9);
            obj.control.axesToolBar   = newPanel(panel,[20 499 570 25,1 1 0 1],[.95 .95 1]);
            obj.control.axesButton(1) = newButtonT(obj.control.axesToolBar,[5 1 25 25,1 0 1 0],obj.GUI.icons.zoomIn,'Zoom in','choosed',true,'callback',@(h,~)obj.cb_axesButton(h,1));
            obj.control.axesButton(2) = newButtonT(obj.control.axesToolBar,[31 1 25 25,1 0 1 0],obj.GUI.icons.zoomOut,'Zoom out','choosed',true,'callback',@(h,~)obj.cb_axesButton(h,2));
            obj.control.axesButton(3) = newButtonT(obj.control.axesToolBar,[57 1 25 25,1 0 1 0],obj.GUI.icons.pan,'Pan','choosed',true,'callback',@(h,~)obj.cb_axesButton(h,3));
            obj.control.axesButton(4) = newButtonT(obj.control.axesToolBar,[83 1 25 25,1 0 1 0],obj.GUI.icons.rotate,'Rotate','choosed',true,'callback',@(h,~)obj.cb_axesButton(h,4));
            newLine(obj.control.axesToolBar,[113 3 1 20,1 0 1 0]);
            x = newButtonT(obj.control.axesToolBar,[118 1 22 25,1 0 1 0],obj.GUI.icons.window,'Show in new figure','callback',@obj.cb_figureMenuButton);
            obj.control.figureMenu = newMenu(panel,[1 1 155 20]);
            obj.control.figureMenu.add([],'Show in new figure','callback',@obj.cb_figureMenu);
            newTip(obj.control.axesToolBar,[1 1 12 25,1 0 1 0],x,'callback',@obj.cb_figureMenuButton);
            obj.control.showButton  = newPopmenu3(obj.control.axesToolBar,[453 1 105 25,0 1 1 0],[{obj.GUI.icons.showPF,obj.GUI.icons.showPS,[],obj.GUI.icons.showTPF,[]},cell(1,length(obj.GUI.metList{1,2})+length(obj.GUI.metList{2,2}))],...
                                                  [{'Pareto Front','Pareto Set','','True PF',''},obj.GUI.metList{2,2}',obj.GUI.metList{1,2}'],'callback',@obj.cb_showbutton);
            obj.control.axes        = newAxes(panel,[82 108 452 360,1 1 1 1]);
            obj.control.playbar     = newPlay(panel,[20 32 570 27,1 1 1 0],'callback',@obj.cb_playbar);
            obj.control.playlabel   = newLabel(panel,[20 7 150 20,1 0 1 0],'','FontSize',9,'HorizontalAlignment','left','ForegroundColor',[.4 .4 .4]);
            obj.control.playStart   = newButtonI(panel,[290 5 30 30,0 0 1 0],obj.GUI.icons.start,'callback',@obj.cb_start);
            obj.control.playPause   = newButtonI(panel,[290 5 30 30,0 0 1 0],obj.GUI.icons.pause,'visible',false,'callback',@obj.cb_pause);
            obj.control.playBack    = newButtonI(panel,[250 5 30 30,0 0 1 0],obj.GUI.icons.back,'callback',@obj.cb_back);
            obj.control.playForward = newButtonI(panel,[330 5 30 30,0 0 1 0],obj.GUI.icons.forward,'callback',@obj.cb_forward);
            obj.control.playStop    = newButtonI(panel,[210 5 30 30,0 0 1 0],obj.GUI.icons.stop,'callback',@obj.cb_stop);
            obj.control.waitlabel   = newLabel(panel,[265 260 100 20,0 0 0 0],'Please wait ... ...','visible',false);
            
            % The fourth panel
            panel = Panel.panels(4);
            newLabel(panel,[10 530 209 15,0 0 0 1],'Result information','ForegroundColor',[.7 .7 .7],'FontSize',9);
            obj.control.resultList = newPopmenu2(panel,[10 490 210 20,0 0 0 1],[.4 .4 1],true,[],'callback',@obj.cb_resultList,'delcallback',@obj.cb_resultListDel);
            obj.control.metricList = newPopmenu2(panel,[10 455 120 20,0 0 0 1],[.1 .1 .1],false,['runtime';obj.GUI.metList{2,2};obj.GUI.metList{1,2}],'callback',@obj.cb_metricList);
            obj.control.metricList.index = 1;
            obj.control.metricLabel = newLabel(panel,[140 453 80 20,0 0 0 1],'','ForegroundColor',[.1 .1 .1],'FontSize',10);
            obj.control.setList     = newLabelM(panel,[10 13 210 427,0 0 1 1],'');
        end
    end
    methods(Static)
        %% Obtain the summary
        function words = summary()
            words = ['Test one multi-objective evolutionary algorithm on a multi-objective ',...
            'optimization problem with specified parameter settings. You can analyse the result and ',...
            'study the performance of the algorithm from several aspects.'];
        end
        %% Callback of the figure in guidance mode
        function guidance(obj)
            obj.GUI.figure.enable        = false;
            obj.control.guideLabel.state = true;
            obj.control.guideIndex       = obj.control.guideIndex + 1;
            switch obj.control.guideIndex
                case 1
                    obj.updateGuideLabel(' Select the algorithm to be executed ',obj.control.setLabel(1));
                    obj.control.algPopmenu.button.state = true;
                case 2
                    obj.updateGuideLabel(' Select the problem to be solved ',obj.control.setLabel(2));
                    obj.control.proPopmenu.button.state = true;
                case 3
                    obj.updateGuideLabel(' Set the parameter values of the algorithm and problem ',obj.control.setPanel);
                case 4
                    obj.updateGuideLabel(' Run the algorithm ',obj.control.runButton);
                case 5
                    obj.updateGuideLabel(' Display the result from different aspects ',obj.control.axesToolBar);
                case 6
                    obj.updateGuideLabel(' Select the result to be displayed ',obj.control.resultList.button);
                otherwise
                    obj.GUI.figure.busy   = false;
                    obj.GUI.figure.enable = true;
                    delete(obj.control.guideLabel);
            end
        end
    end
    methods(Access = private)
        %% Update the parameter setting list
        function cb_apoPopmenu(obj,type)
            switch type
                case 'alg'
                    obj.control.setPanel.update('replace',obj.control.algPopmenu.string,[.2 .4 .7],1,1)
                case 'pro'
                    obj.control.setPanel.update('replace',obj.control.proPopmenu.string,[.9 .5 .2],-1,2);
            end
        end
        %% Run the algorithm
        function cb_run(obj,hObject,eventdata)
            % Read all the parameters
            items     = obj.control.setPanel.items;
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
            % Generate the input for executing algorithm
            len    = cumsum(arrayfun(@(S)length(S.edits),items));
            Inputs = {'-N',Parameter{len(1)+1},'-M',Parameter{len(1)+2},'-D',Parameter{len(1)+3},'-evaluation',Parameter{len(1)+4},...
                      '-algorithm',[{str2func(items(1).name)},Parameter(1:len(1))'],...
                      '-problem',[{str2func(items(2).name)},Parameter(len(1)+5:len(2))']};
            % Update the data
            try
                obj.data = [obj.data,module_test_result(sprintf('%s on %s',items(1).name,items(2).name),...
                                                        GLOBAL(Inputs{:},'-outputFcn',@obj.outputFcn))];
                obj.control.resultList.add(obj.data(end).name);
                obj.control.resultList.index = length(obj.control.resultList.items);
            catch err
                errordlg(err.message,get(gcf,'Name'),'modal'); beep;
                return;
            end
            % Update the state of controls
            obj.GUI.menu.enable           = false;
            obj.panel.panels(1).enable    = false;
            obj.panel.panels(2).enable    = false;
            obj.panel.panels(4).enable    = false;
            obj.control.playStart.visible = false;
            obj.control.playPause.visible = true;
            obj.control.playStop.value    = false;
            [obj.control.showButton.menu.items(3:end).enable] = deal(false);
            obj.control.showButton.index  = 1;
            % Call the algorithm
            try
                obj.data(end).Global.Start();
                obj.cb_stop();
            catch err
                obj.cb_stop();
                rethrow(err);
            end
        end
        %% Axis controlling
        function cb_axesButton(obj,hObject,eventdata)
            [obj.control.axesButton([1:eventdata-1,eventdata+1:end]).value] = deal(false);
            obj.control.axes.mouseKind = hObject.value*eventdata;
        end
        %% Show the figure menu
        function cb_figureMenuButton(obj,hObject,eventdata)
            [obj.control.figureMenu.items.enable] = deal(~isempty(obj.control.axes.handle.Children));
            obj.control.figureMenu.show();
        end
        %% Show the axis in new figure
        function cb_figureMenu(obj,hObject,eventdata)
            persistent nFigure;
            if isempty(nFigure)
                nFigure = 1;
            else
                nFigure = nFigure + 1;
            end
            % Create a new figure
            newFigure = figure('NumberTitle','off','Name',sprintf('Figure %d',nFigure),'DeleteFcn',@(obj,~)obj.UserData.parent.del([obj.UserData.parent.items.handle]==obj.UserData.handle));
            % Create a new menu item
            obj.control.figureMenu.add([],['Show in ',newFigure.Name],'callback',@obj.cb_figureMenuItem,'UserData',newFigure);
            % Associate the new figure with the new menu item
            newFigure.UserData = obj.control.figureMenu.items(end);
            % Copy the axis to the new figure
            newAxes = copyobj(obj.control.axes.handle,newFigure);
            set(newAxes,'Units','normalized','Position',[0.1300 0.1100 0.7750 0.8150]);
            str = newAxes.Title.String;
            newAxes.Title.String = str(strfind(str,' on ')+4:end);
            legend(str(1:strfind(str,' on ')-1),'Location','NorthEast');
            newAxes.Legend.Interpreter = 'none';
        end
        %% Show the axis in existing figure
        function cb_figureMenuItem(obj,hObject,eventdata)
            figure(hObject.handle.UserData);
            currentAxes = gca;
            color = rand(1,3);
            set(copyobj(obj.control.axes.handle.Children,currentAxes),'Color',color,'Markerfacecolor',color);
            str = obj.control.axes.handle.Title.String;
            currentAxes.Legend.String{end} = str(1:strfind(str,' on ')-1);
            axis(gca,'tight');
        end
        %% Running - start
        function cb_start(obj,hObject,eventdata)
            obj.control.playPause.moved = true;
            if obj.panel.panels(1).enable
                obj.cb_run();
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
        %% Running - back
        function cb_back(obj,hObject,eventdata)
            if obj.control.resultList.index > 0
                if obj.control.playPause.visible
                    obj.control.playStart.visible = true;
                    obj.control.playPause.visible = false;
                end
                obj.control.playbar.value = obj.control.playbar.value-obj.control.playbar.maxvalue./(size(obj.data(obj.control.resultList.index).result,1)-1);
            end
        end
        %% Running - forward
        function cb_forward(obj,hObject,eventdata)
            if obj.control.resultList.index > 0
                if obj.control.playPause.visible
                    obj.control.playStart.visible = true;
                    obj.control.playPause.visible = false;
                end
                if obj.control.playbar.value < obj.control.playbar.maxvalue || obj.panel.panels(1).enable
                    obj.control.playbar.value = obj.control.playbar.value+obj.control.playbar.maxvalue./(size(obj.data(obj.control.resultList.index).result,1)-1);
                else
                    obj.control.playStart.handle.UserData = rand;
                end
            end
        end
        %% Running - stop
        function cb_stop(obj,hObject,eventdata)
            if ~obj.panel.panels(1).enable
                % Update the state of controls
                obj.control.playStart.handle.UserData = rand;
                obj.GUI.menu.enable           = true;
                obj.panel.panels(1).enable    = true;
                obj.panel.panels(2).enable    = true;
                obj.panel.panels(4).enable    = true;
                obj.control.playStart.visible = true;
                obj.control.playPause.visible = false;
                [obj.control.showButton.menu.items(3:end).enable] = deal(true);
                % Show the last result
                obj.data(end).finalset();
                obj.data(end).playmaxvalue   = obj.control.playbar.maxvalue;
                obj.control.resultList.index = length(obj.control.resultList.items);
            end
        end
        %% Show the population specified by the play bar
        function cb_playbar(obj,hObject,eventdata)
            current = obj.control.resultList.index;
            if current > 0 && ~isempty(obj.data(current).result)
                index = 1 + round(obj.control.playbar.value./obj.control.playbar.maxvalue*(size(obj.data(current).result,1)-1));
                obj.control.playlabel.handle.String = sprintf('%d/%d',obj.data(current).result{index,1},obj.data(current).result{end,1});
                switch obj.control.showButton.handle.String
                    case 'Pareto Front'
                        obj.control.axes.draw(obj.data(current).result{index,2});
                    case 'Pareto Set'
                        obj.control.axes.draw(obj.data(current).result{index,3});
                    case 'True PF'
                        if ~obj.control.axes.viewFix
                            obj.control.axes.draw(obj.data(current).Global.PF);
                        end
                    otherwise
                        if ~obj.control.axes.viewFix
                            metric = obj.control.showButton.handle.String;
                            if ~isfield(obj.data(current).metricline,metric)
                                obj.GUI.figure.busy = true;
                                obj.control.waitlabel.visible = true;
                                drawnow();
                                selected = unique(ceil(linspace(1,size(obj.data(current).result,1),10)));
                                obj.data(current).metricline.(metric)(:,1) = cell2mat(obj.data(current).result(selected,1));
                                obj.data(current).metricline.(metric)(:,2) = cellfun(@(S)obj.Metric(str2func(metric),S,obj.data(current).Global.PF),obj.data(current).result(selected,2));
                                obj.GUI.figure.busy = false;
                                obj.control.waitlabel.visible = false;
                            end
                            obj.control.axes.draw(obj.data(current).metricline.(metric),'-k.','LineWidth',1.5,'MarkerSize',10);
                            xlabel('Number of evaluations');
                            ylabel(metric);
                        end
                end
                if ~isempty(obj.control.axes.handle.Children)
                    obj.control.axes.viewFix = true;
                end
            else
                obj.control.playlabel.handle.String = '';
                cla(obj.control.axes.handle);
            end
        end
        %% Show the specified type of result on the axis
        function cb_showbutton(obj,hObject,eventdata)
            obj.control.axes.viewFix = false;
            obj.cb_playbar();
        end
        %% Show the specified recorded result
        function cb_resultList(obj,hObject,eventdata)
            current = obj.control.resultList.index;
            if current > 0
                obj.control.playbar.maxvalue      = obj.data(current).playmaxvalue;
                obj.control.setList.handle.String = obj.data(current).setting;
                obj.cb_metricList();
                obj.cb_showbutton();
                title(obj.data(current).name,'Interpreter','none');
            else
                obj.control.playbar.maxvalue      = 0;
                obj.control.setList.handle.String = '';
                obj.cb_metricList();
                obj.cb_showbutton();
                title('');
            end
        end
        %% Update the result list after deleting one result
        function cb_resultListDel(obj,hObject,eventdata)
            obj.data(eventdata) = [];
        end
        %% Calculate the metric value of the final population
        function cb_metricList(obj,hObject,eventdata)
            current = obj.control.resultList.index;
            if current > 0 && ~isempty(obj.data(current).result)
                metric = obj.control.metricList.string;
                if ~isfield(obj.data(current).metric,metric)
                    obj.data(current).metric.(metric) = obj.Metric(str2func(metric),obj.data(current).result{end,2},obj.data(current).Global.PF);
                end
                obj.control.metricLabel.handle.String = sprintf('%.4e',obj.data(current).metric.(metric));
            else
                obj.control.metricLabel.handle.String = '';
            end
        end
        %% Output function
        function outputFcn(obj,Global)
            Population = Global.result{end};
            if ~isempty(Population)
                Feasible = all(Population.cons<=0,2);
                obj.data(end).result = [obj.data(end).result;{Global.evaluated},{Population(Feasible).objs},{Population(Feasible).decs}];
            else
                obj.data(end).result = [obj.data(end).result;{Global.evaluated},{[]},{[]}];
            end
            obj.control.playbar.maxvalue = size(obj.data(end).result,1)./(size(obj.data(end).result,1)+max(0,Global.maxgen-Global.gen));
            obj.control.playbar.value    = obj.control.playbar.maxvalue;
            if obj.control.playStop.value
                error('GLOBAL:Termination','Algorithm has terminated');
            end
            if obj.control.playStart.visible
            	waitfor(obj.control.playStart.handle,'UserData');
            end
            if obj.control.playStop.value
                error('GLOBAL:Termination','Algorithm has terminated');
            end
        end
    end
end