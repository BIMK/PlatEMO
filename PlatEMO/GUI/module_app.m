classdef module_app < handle
%module_test - Application module.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        GUI;                % The GUI object
        app  = struct(); 	% All the components
        pro  = [];          % Components of problem functions
        con  = [];          % Components of constraint functions
        data = {};          % All the results
    end
    methods(Access = ?GUI)
        %% Constructor
        function obj = module_app(GUI)
            % The main grid
            obj.GUI = GUI;
            obj.app.maingrid = GUI.APP(3,1,uigridlayout(obj.GUI.app.maingrid,'RowHeight',{20,30,'1x'},'ColumnWidth',{'2x',5,1,'1.2x',1,'3x'},'Padding',[5 5 5 5],'RowSpacing',5,'ColumnSpacing',0,'BackgroundColor','w'));
            obj.app.label(1) = GUI.APP(1,1,uilabel(obj.app.maingrid,'Text','Problem definition','HorizontalAlignment','center','FontSize',11,'FontColor',[.5 .5 .5]));
            obj.app.label(2) = GUI.APP(1,4,uilabel(obj.app.maingrid,'Text','Algorithm selection','HorizontalAlignment','center','FontSize',11,'FontColor',[.5 .5 .5]));
            obj.app.label(3) = GUI.APP(1,6,uilabel(obj.app.maingrid,'Text','Result display','HorizontalAlignment','center','FontSize',11,'FontColor',[.5 .5 .5]));
            GUI.APP([1 3],3,uipanel(obj.app.maingrid,'BackgroundColor',[.8 .8 .8]));
            GUI.APP([1 3],5,uipanel(obj.app.maingrid,'BackgroundColor',[.8 .8 .8]));
            
            % The first panel
            obj.app.grid(1)    = GUI.APP(2,1,uigridlayout(obj.app.maingrid,'RowHeight',{'1x'},'ColumnWidth',{'1x','1x','1.1x','1.1x','1x'},'Padding',[5 5 5 5],'RowSpacing',5,'ColumnSpacing',5,'BackgroundColor',[.95 .95 1]));
            obj.app.buttonA(1) = GUI.APP(1,1,uibutton(obj.app.grid(1),'Text','+ objective','BackgroundColor','w','Tooltip','Increase the number of objectives','Interruptible','off','BusyAction','cancel','ButtonpushedFcn',{@obj.cb_updateProblem,1}));
            obj.app.buttonA(2) = GUI.APP(1,2,uibutton(obj.app.grid(1),'Text','- objective','BackgroundColor','w','Tooltip','Decrease the number of objectives','Interruptible','off','BusyAction','cancel','Enable','off','ButtonpushedFcn',{@obj.cb_updateProblem,2}));
            obj.app.buttonA(3) = GUI.APP(1,3,uibutton(obj.app.grid(1),'Text','+ constraint','BackgroundColor','w','Tooltip','Increase the number of constraints','Interruptible','off','BusyAction','cancel','ButtonpushedFcn',{@obj.cb_updateProblem,3}));
            obj.app.buttonA(4) = GUI.APP(1,4,uibutton(obj.app.grid(1),'Text','- constraint','BackgroundColor','w','Tooltip','Decrease the number of constraints','Interruptible','off','BusyAction','cancel','ButtonpushedFcn',{@obj.cb_updateProblem,4}));
            obj.app.buttonA(5) = GUI.APP(1,5,uibutton(obj.app.grid(1),'Text','Validation','BackgroundColor',[.95 .95 1],'Tooltip','Check the validity of the problem','Interruptible','off','BusyAction','cancel','ButtonpushedFcn',@obj.cb_validation));
            
            % The second panel
            obj.app.grid(2)    = GUI.APP(3,1,uigridlayout(obj.app.maingrid,'RowHeight',repmat({22},1,15),'ColumnWidth',{75,'1x',25,20},'Padding',[5 10 5 5],'RowSpacing',5,'ColumnSpacing',5,'Scrollable','on','BackgroundColor','w'));
            obj.app.titleB(1)  = GUI.APP(1,[1 4],uilabel(obj.app.grid(2),'Text','Decision vector','FontSize',15,'FontColor',[.9 .5 .2],'FontWeight','bold'));
            tempGrid           = GUI.APP(2,[1 4],uigridlayout(obj.app.grid(2),'RowHeight',{'1x'},'ColumnWidth',{75,40,100,60,'1x'},'Padding',[0 0 0 0],'RowSpacing',0,'ColumnSpacing',5,'BackgroundColor','w'));
            obj.app.labelB(1)  = GUI.APP(1,1,uilabel(tempGrid,'Text','x =','HorizontalAlignment','right'));
            obj.app.editB(1)   = GUI.APP(1,2,uieditfield(tempGrid,'numeric','Value',10,'limits',[1 inf],'RoundFractionalValues','on','ValueChangedFcn',@obj.cb_updateFilter));
            obj.app.dropB      = GUI.APP(1,3,uidropdown(tempGrid,'BackgroundColor','w','Items',{'real','binary','permutation'},'ItemsData',1:3,'Value',1,'ValueChangedFcn',@obj.cb_updateEncoding));
            obj.app.labelB(2)  = GUI.APP(1,4,uilabel(tempGrid,'Text','variables'));
            obj.app.titleB(2)  = GUI.APP(3,[1 4],uilabel(obj.app.grid(2),'Text','Decision space','FontSize',15,'FontColor',[.9 .5 .2],'FontWeight','bold'));
            obj.app.labelB(3)  = GUI.APP(4,1,uilabel(obj.app.grid(2),'Text','x >=','HorizontalAlignment','right'));
            obj.app.editB(2)   = GUI.APP(4,[2 3],uieditfield(obj.app.grid(2),'Value','zeros(1,10)-1','Tooltip','Lower bound of each decision variable'));
            obj.app.labelB(4)  = GUI.APP(5,1,uilabel(obj.app.grid(2),'Text','x <=','HorizontalAlignment','right'));
            obj.app.editB(3)   = GUI.APP(5,[2 3],uieditfield(obj.app.grid(2),'Value','zeros(1,10)+1','Tooltip','Upper bound of each decision variable'));
            obj.app.titleB(3)  = GUI.APP(6,[1 4],uilabel(obj.app.grid(2),'Text','Dataset','FontSize',15,'FontColor',[.9 .5 .2],'FontWeight','bold'));
            obj.app.labelB(5)  = GUI.APP(7,1,uilabel(obj.app.grid(2),'Text','data =','HorizontalAlignment','right'));
            obj.app.editB(4)   = GUI.APP(7,[2 3],uieditfield(obj.app.grid(2),'Value','eye(10)','Tooltip','Dataset of the problem'));
            obj.app.buttonB(1) = GUI.APP(7,4,uibutton(obj.app.grid(2),'Text','...','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_loadFunction,obj.app.editB(4),'*.mat'}));
            obj.app.titleB(4)  = GUI.APP(8,[1 4],uilabel(obj.app.grid(2),'Text','Initialization function','FontSize',15,'FontColor',[.9 .5 .2],'FontWeight','bold'));
            obj.app.labelB(6)  = GUI.APP(9,1,uilabel(obj.app.grid(2),'Text','I(n,data) =','HorizontalAlignment','right'));
            obj.app.editB(5)   = GUI.APP(9,[2 3],uieditfield(obj.app.grid(2),'Value','rand(n,10)*2-1','Tooltip','Function for generating an initial population with size n'));
            obj.app.buttonB(2) = GUI.APP(9,4,uibutton(obj.app.grid(2),'Text','...','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_loadFunction,obj.app.editB(5),'*.m'}));
            obj.app.titleB(5)  = GUI.APP(10,[1 4],uilabel(obj.app.grid(2),'Text','Repair function','FontSize',15,'FontColor',[.9 .5 .2],'FontWeight','bold'));
            obj.app.labelB(7)  = GUI.APP(11,1,uilabel(obj.app.grid(2),'Text','R(x,data) =','HorizontalAlignment','right'));
            obj.app.editB(6)   = GUI.APP(11,[2 3],uieditfield(obj.app.grid(2),'Value','max(min(x,1),-1)','Tooltip','Function for repairing invalid solution'));
            obj.app.buttonB(3) = GUI.APP(11,4,uibutton(obj.app.grid(2),'Text','...','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_loadFunction,obj.app.editB(6),'*.m'}));
            obj.app.titleB(6)  = GUI.APP(12,[1 4],uilabel(obj.app.grid(2),'Text','Objective functions','FontSize',15,'FontColor',[.9 .5 .2],'FontWeight','bold'));
            obj.app.titleB(7)  = GUI.APP(14,[1 4],uilabel(obj.app.grid(2),'Text','Constraint functions','FontSize',15,'FontColor',[.9 .5 .2],'FontWeight','bold'));
            obj.pro(1).label   = GUI.APP(13,1,uilabel(obj.app.grid(2),'Text','f1(x,data) =','HorizontalAlignment','right'));
            obj.pro(1).edit    = GUI.APP(13,[2 3],uieditfield(obj.app.grid(2),'Value','mean(x*data)','Tooltip','Objective function to be minimized'));
            obj.pro(1).button  = GUI.APP(13,4,uibutton(obj.app.grid(2),'Text','...','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_loadFunction,obj.pro(1).edit,'*.m'}));
            obj.con(1).labelA  = GUI.APP(15,1,uilabel(obj.app.grid(2),'Text','g1(x,data) =','HorizontalAlignment','right'));
            obj.con(1).edit    = GUI.APP(15,2,uieditfield(obj.app.grid(2),'Value','mean(x)','Tooltip','Constraint function to be satisfied'));
            obj.con(1).labelB  = GUI.APP(15,3,uilabel(obj.app.grid(2),'Text','<= 0'));
            obj.con(1).button  = GUI.APP(15,4,uibutton(obj.app.grid(2),'Text','...','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_loadFunction,obj.con(1).edit,'*.m'}));
            
            % The third panel
            obj.app.grid(3)    = GUI.APP([2 3],4,uigridlayout(obj.app.maingrid,'RowHeight',{16,20,16,20,16,20,20,25,'1.5x','1x'},'ColumnWidth',{'1x','1x','1x'},'Padding',[8 10 8 0],'RowSpacing',5,'ColumnSpacing',5,'BackgroundColor','w'));
            obj.app.labelC(1)  = GUI.APP(1,[1 3],uilabel(obj.app.grid(3),'Text','Number of objectives','VerticalAlignment','bottom','FontSize',12,'FontColor',[.15 .6 .2],'FontWeight','bold'));
            obj.app.stateC(1)  = GUI.APP(2,1,uibutton(obj.app.grid(3),'state','Text','single','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','The problem has a single objective','ValueChangedFcn',{@obj.cb_filter,1}));
            obj.app.stateC(2)  = GUI.APP(2,2,uibutton(obj.app.grid(3),'state','Text','multi','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','The problem has 2 or 3 objectives','ValueChangedFcn',{@obj.cb_filter,2}));
            obj.app.stateC(3)  = GUI.APP(2,3,uibutton(obj.app.grid(3),'state','Text','many','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','The problem has more than 3 objectives','ValueChangedFcn',{@obj.cb_filter,3}));
            obj.app.labelC(2)  = GUI.APP(3,[1 3],uilabel(obj.app.grid(3),'Text','Encoding scheme','VerticalAlignment','bottom','FontSize',12,'FontColor',[.15 .6 .2],'FontWeight','bold'));
            obj.app.stateC(4)  = GUI.APP(4,1,uibutton(obj.app.grid(3),'state','Text','real','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','The decision variables are real values','ValueChangedFcn',{@obj.cb_filter,4}));
            obj.app.stateC(5)  = GUI.APP(4,2,uibutton(obj.app.grid(3),'state','Text','binary','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','The decision variables are binary values','ValueChangedFcn',{@obj.cb_filter,5}));
            obj.app.stateC(6)  = GUI.APP(4,3,uibutton(obj.app.grid(3),'state','Text','permutation','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','The decision vector is a permutation','ValueChangedFcn',{@obj.cb_filter,6}));
            obj.app.labelC(3)  = GUI.APP(5,[1 3],uilabel(obj.app.grid(3),'Text','Special difficulties','VerticalAlignment','bottom','FontSize',12,'FontColor',[.15 .6 .2],'FontWeight','bold'));
            obj.app.stateC(7)  = GUI.APP(6,1,uibutton(obj.app.grid(3),'state','Text','large','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','The problem has more than 100 decision variables','ValueChangedFcn',{@obj.cb_filter,7}));
            obj.app.stateC(8)  = GUI.APP(6,2,uibutton(obj.app.grid(3),'state','Text','constrained','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','The problem has constraints','ValueChangedFcn',{@obj.cb_filter,8}));
            obj.app.stateC(9)  = GUI.APP(6,3,uibutton(obj.app.grid(3),'state','Text','expensive','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','The objectives are computationally time-consuming','ValueChangedFcn',{@obj.cb_filter,9}));
            obj.app.stateC(10) = GUI.APP(7,1,uibutton(obj.app.grid(3),'state','Text','multimodal','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','The objectives are multimodal','ValueChangedFcn',{@obj.cb_filter,10}));
            obj.app.stateC(11) = GUI.APP(7,2,uibutton(obj.app.grid(3),'state','Text','sparse','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','Most decision variables of the optimal solutions are zero','ValueChangedFcn',{@obj.cb_filter,11}));
            obj.app.stateC(12) = GUI.APP(7,3,uibutton(obj.app.grid(3),'state','Text','preference','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','Only searching for preferred regions of the Pareto front','ValueChangedFcn',{@obj.cb_filter,12}));
            obj.app.labelC(4)  = GUI.APP(8,[1 2],uilabel(obj.app.grid(3),'Text','Algorithms','VerticalAlignment','bottom','FontSize',13,'FontColor',[.2 .4 .7],'FontWeight','bold'));
            obj.app.labelC(5)  = GUI.APP(8,3,uilabel(obj.app.grid(3),'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',10,'FontColor',[.2 .4 .7]));
            obj.app.listC      = GUI.APP(9,[1 3],uilistbox(obj.app.grid(3),'FontColor',[.2 .4 .7],'ValueChangedFcn',{@obj.cb_updateList,3}));
            obj.app.listD      = uilist(obj.app.grid(3),obj.GUI.app.figure,obj.GUI.iconFolder);
            obj.app.tipD       = GUI.APP(10,[1 3],uilabel(obj.app.grid(3),'Text','Select an algorithm','HorizontalAlignment','center'));
            obj.app.listD.grid.Padding = [0,0,0,0];
            GUI.APP(10,[1 3],obj.app.listD.grid);
            
            % The fourth panel
            obj.app.grid(4) = GUI.APP([2 3],6,uigridlayout(obj.app.maingrid,'RowHeight',{'1x',40,30},'ColumnWidth',{20,150,'1x','1x',150,20},'Padding',[15 10 15 0],'RowSpacing',5,'BackgroundColor','w'));
            obj.app.axes    = GUI.APP(1,[2 5],uiaxes(obj.app.grid(4),'BackgroundColor','w','Box','on'));
            obj.app.waittip = GUI.APP(1,[2 5],uilabel(obj.app.grid(4),'HorizontalAlignment','center','Text','                 Please wait ... ...','Visible','off'));
            obj.app.menu(1) = uicontext2(obj.GUI.app.figure,@obj.cb_slider);
            obj.app.menu(1).add('Population (obj.)',false);
            obj.app.menu(1).add('Population (dec.)',true);
            obj.app.menu(1).add('HV',false);
            obj.app.menu(1).add('Feasible rate',false);
            obj.app.menu(1).flush();
            obj.app.menu(2) = uicontext2(obj.GUI.app.figure,@obj.cb_slider);
            obj.app.menu(2).add('Population (dec.)',true);
            obj.app.menu(2).add('Minimum value',false);
            obj.app.menu(2).add('Feasible rate',false);
            obj.app.menu(2).flush();
            tempTb = axtoolbar(obj.app.axes(1),{'rotate','pan','zoomin','zoomout'});
            obj.app.toolD(1)   = axtoolbarbtn(tempTb,'push','Icon',fullfile(obj.GUI.iconFolder,'gif.png'),'Tooltip','Save the evolutionary process to gif','ButtonPushedFcn',@obj.cb_toolbutton1);
            obj.app.toolD(2)   = axtoolbarbtn(tempTb,'push','Icon',fullfile(obj.GUI.iconFolder,'newfigure.png'),'Tooltip','Open in new figure and save to workspace','ButtonPushedFcn',@obj.cb_toolbutton2);
            obj.app.toolD(3)   = axtoolbarbtn(tempTb,'push','Icon',fullfile(obj.GUI.iconFolder,'datasource.png'),'Tooltip','Data source','ButtonPushedFcn',@obj.cb_toolbutton3);
            obj.app.slider     = GUI.APP(2,[1 6],uislider(obj.app.grid(4),'Limits',[0 1],'MajorTicks',0:0.25:1,'MajorTickLabels',{'0%','25%','50%','75%','100%'},'MinorTicks',0:0.01:1,'ValueChangedFcn',@obj.cb_slider));
            obj.app.labelD     = GUI.APP(3,[5 6],uilabel(obj.app.grid(4),'Text','','HorizontalAlignment','right'));
            obj.app.buttonD(1) = GUI.APP(3,3,uibutton(obj.app.grid(4),'push','Text','Start','FontSize',16,'Enable','off','ButtonpushedFcn',@obj.cb_start));
            obj.app.buttonD(2) = GUI.APP(3,4,uibutton(obj.app.grid(4),'push','Text','Stop','FontSize',16,'Enable','off','ButtonpushedFcn',@obj.cb_stop));
            
            % Initialization
            obj.cb_updateFilter();
        end
    end
    methods(Access = private)
        %% Update the algorithms in the lists
        function cb_filter(obj,~,~,index)
            if nargin > 3
                if index < 4
                    [obj.app.stateC(1:3).Value] = deal(0);
                    obj.app.stateC(index).Value = 1;
                elseif index < 7
                    [obj.app.stateC(4:6).Value] = deal(0);
                    obj.app.stateC(index).Value = 1;
                end
            end
            filter = [obj.app.stateC.Value];
            func   = @(s)all(any(repmat([true,filter],size(s,1),1)&s,2)) && all((any(s(:,2:end),1)&filter)==filter);
            show   = cellfun(func,obj.GUI.algList(:,1));
            obj.app.listC.Items = ['(Open File)';obj.GUI.algList(show,2)];
            obj.app.listC.Value = {};
            obj.app.labelC(5).Text = sprintf('%d / %d',sum(show),length(show));
            obj.app.listD.del(1);
            obj.app.tipD.Visible      = 'on';
            obj.app.buttonD(1).Enable = 'off';
        end
        %% Update the parameter list
        function cb_updateList(obj,~,~,type)
            filename = obj.app.listC.Value;
            if contains(filename,'Open File')
                [file,path] = uigetfile('*.m','');
                if file ~= 0
                    try
                        filename = fullfile(path,file);
                        f   = fopen(filename);
                        str = fgetl(f);
                        fclose(f);
                        assert(contains(str,'< ALGORITHM'));
                    catch
                        uialert(obj.GUI.app.figure,'The selected file is not a subclass of ALGORITHM.','Error');
                        return;
                    end
                else
                    return;
                end
            else
                filename = [filename,'.m'];
            end
            obj.app.tipD.Visible      = 'off';
            obj.app.buttonD(1).Enable = 'on';
            obj.app.listD.del(1);
            obj.app.listD.add(filename,type);
            obj.app.listD.flush();
        end
        %% Update the objective and constraint functions
        function cb_updateProblem(obj,~,~,index)
            k = 0;
            switch index
                case 1
                    item.label  = GUI.APP(length(obj.pro)+13,1,uilabel(obj.app.grid(2),'Text',sprintf('f%d(x,data) =',length(obj.pro)+1),'HorizontalAlignment','right'));
                    item.edit   = GUI.APP(length(obj.pro)+13,[2 3],uieditfield(obj.app.grid(2),'Value','mean(x*data)','Tooltip','Objective function to be minimized'));
                    item.button = GUI.APP(length(obj.pro)+13,4,uibutton(obj.app.grid(2),'Text','...','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_loadFunction,item.edit,'*.m'}));
                    obj.pro     = [obj.pro,item];
                    k = 1;
                case 2
                    if length(obj.pro) > 1
                        delete(obj.pro(end).label);
                        delete(obj.pro(end).edit);
                        delete(obj.pro(end).button);
                        obj.pro(end) = [];
                        k = -1;
                    end
                case 3
                    item.labelA = GUI.APP(length(obj.pro)+14+length(obj.con),1,uilabel(obj.app.grid(2),'Text',sprintf('g%d(x,data) =',length(obj.con)+1),'HorizontalAlignment','right'));
                    item.edit   = GUI.APP(length(obj.pro)+14+length(obj.con),2,uieditfield(obj.app.grid(2),'Value','mean(x)','Tooltip','Constraint function to be satisfied'));
                    item.labelB = GUI.APP(length(obj.pro)+14+length(obj.con),3,uilabel(obj.app.grid(2),'Text','<= 0'));
                    item.button = GUI.APP(length(obj.pro)+14+length(obj.con),4,uibutton(obj.app.grid(2),'Text','...','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_loadFunction,item.edit,'*.m'}));
                    obj.con     = [obj.con,item];
                case 4
                    if ~isempty(obj.con)
                        delete(obj.con(end).labelA);
                        delete(obj.con(end).edit);
                        delete(obj.con(end).labelB);
                        delete(obj.con(end).button);
                        obj.con(end) = [];
                    end
            end
            obj.app.buttonA(2).Enable = length(obj.pro) > 1;
            obj.app.buttonA(4).Enable = ~isempty(obj.con);
            if k ~= 0
                obj.app.titleB(7).Layout.Row = obj.app.titleB(7).Layout.Row + k;
                for i = 1 : length(obj.con)
                    obj.con(i).labelA.Layout.Row = obj.con(i).labelA.Layout.Row + k;
                    obj.con(i).edit.Layout.Row   = obj.con(i).edit.Layout.Row + k;
                    obj.con(i).labelB.Layout.Row = obj.con(i).labelB.Layout.Row + k;
                    obj.con(i).button.Layout.Row = obj.con(i).button.Layout.Row + k;
                end
            end
            obj.app.grid(2).RowHeight = repmat({22},1,13+length(obj.pro)+length(obj.con));
            if isempty(obj.con)
                obj.app.titleB(7).Text = 'No constraint';
            else
                obj.app.titleB(7).Text = 'Constraint functions';
            end
            obj.cb_updateFilter();
        end
        %% Update the initialization and repair functions
        function cb_updateEncoding(obj,~,~)
            set(obj.app.editB(2:3),'Enable',obj.app.dropB.Value==1);
            if obj.app.dropB.Value == 1
                obj.app.editB(2).Value = sprintf('zeros(1,%d)-1',obj.app.editB(1).Value);
                obj.app.editB(3).Value = sprintf('zeros(1,%d)+1',obj.app.editB(1).Value);
            elseif obj.app.dropB.Value == 2
                obj.app.editB(2).Value = sprintf('zeros(1,%d)',obj.app.editB(1).Value);
                obj.app.editB(3).Value = sprintf('ones(1,%d)',obj.app.editB(1).Value);
            elseif obj.app.dropB.Value == 3
                obj.app.editB(2).Value = sprintf('zeros(1,%d)+1',obj.app.editB(1).Value);
                obj.app.editB(3).Value = sprintf('zeros(1,%d)+%d',obj.app.editB(1).Value,obj.app.editB(1).Value);
            end
            str = {sprintf('rand(n,%d)*2-1',obj.app.editB(1).Value),sprintf('rand(n,%d)<0.5',obj.app.editB(1).Value),sprintf('feval2(@sort,rand(n,%d),2)',obj.app.editB(1).Value)};
            obj.app.editB(5).Value = str{obj.app.dropB.Value};
            str = {'max(min(x,1),-1)','round(x)','feval2(@sort,feval2(@sort,x))'};
            obj.app.editB(6).Value = str{obj.app.dropB.Value};
            obj.cb_updateFilter();
        end
        %% Validate the problem
        function PRO = cb_validation(obj,para,~)
            msg = ', please refer to the command window for details.';
            try
                lower = str2num(obj.app.editB(2).Value);
                assert(~isempty(lower),'The lower bound cannot be empty');
            catch err
                uialert(obj.GUI.app.figure,['The lower bound of decision space is illegal',msg],'Invalid problem');
                rethrow(err);
            end
            try
                upper = str2num(obj.app.editB(3).Value);
                assert(~isempty(upper),'The upper bound cannot be empty');
            catch err
                uialert(obj.GUI.app.figure,['The upper bound of decision space is illegal',msg],'Invalid problem');
                rethrow(err);
            end
            try
                if ~isempty(regexp(obj.app.editB(4).Value,'^<.+>$','once'))
                    dataset = load(obj.app.editB(4).Value(2:end-1),'-mat');
                else
                    dataset = str2num(obj.app.editB(4).Value);
                end
            catch err
                uialert(obj.GUI.app.figure,['The dataset is illegal',msg],'Invalid problem');
                rethrow(err);
            end
            try
                initFcn = GetFcn(obj.app.editB(5).Value,0);
                PopDec  = initFcn(1,dataset);
            catch err
                uialert(obj.GUI.app.figure,['The initialization function is illegal',msg],'Invalid problem');
                rethrow(err);
            end
            try
                decFcn = GetFcn(obj.app.editB(6).Value);
                PopDec = decFcn(PopDec,dataset);
            catch err
                uialert(obj.GUI.app.figure,['The repair function is illegal',msg],'Invalid problem');
                rethrow(err);
            end
            try
                objFcn = cell(1,length(obj.pro));
                for i = 1 : length(obj.pro)
                    objFcn{i} = GetFcn(obj.pro(i).edit.Value);
                    objFcn{i}(PopDec,dataset);
                end
            catch err
                uialert(obj.GUI.app.figure,['Objective function f',num2str(i),' is illegal',msg],'Invalid problem');
                rethrow(err);
            end
            try
                conFcn = cell(1,length(obj.con));
                for i = 1 : length(obj.con)
                    conFcn{i} = GetFcn(obj.con(i).edit.Value);
                    conFcn{i}(PopDec,dataset);
                end
            catch err
                uialert(obj.GUI.app.figure,['Constraint function g',num2str(i),' is illegal',msg],'Invalid problem');
                rethrow(err);
            end
            if nargout == 0
                uialert(obj.GUI.app.figure,'The problem definition is valid.','Valid problem','Icon','success');
            else
                PRO = UserProblem('N',para{1},'maxFE',para{2},'encoding',obj.app.dropB.Items{obj.app.dropB.Value},'lower',lower,'upper',upper,'initFcn',initFcn,'decFcn',decFcn,'objFcn',objFcn,'conFcn',conFcn,'parameter',dataset);
            end
        end
        %% Load dataset or function
        function cb_loadFunction(obj,~,~,h,ext)
            if length(h.Value)>2 && h.Value(1)=='<' && h.Value(end)=='>'
                [file,folder] = uigetfile(ext,'',fileparts(h.Value(2:end-1)));
            else
                [file,folder] = uigetfile(ext,'',cd);
            end
            if file ~= 0
                h.Value = sprintf('<%s>',fullfile(folder,file));
            end
        end
        %% Update the filter
        function cb_updateFilter(obj,~,~)
            oldState = [obj.app.stateC.Value];
            obj.app.stateC(1).Value = length(obj.pro) < 2;
            obj.app.stateC(2).Value = length(obj.pro) > 1 && length(obj.pro) < 4;
            obj.app.stateC(3).Value = length(obj.pro) > 3;
            [obj.app.stateC(4:6).Value] = deal(0);
            obj.app.stateC(obj.app.dropB.Value+3).Value = 1;
            obj.app.stateC(7).Value = obj.app.dropB.Value~=2 && obj.app.editB(1).Value>99;
            obj.app.stateC(8).Value = ~isempty(obj.con);
            if any(oldState~=[obj.app.stateC.Value])
                obj.cb_filter();
            end
        end
        %% Start the execution
        function cb_start(obj,~,~)
            if strcmp(obj.app.buttonD(1).Text,'Pause')
                obj.app.buttonD(1).Text = 'Continue';
            elseif strcmp(obj.app.buttonD(1).Text,'Continue')
                obj.app.buttonD(1).Text = 'Pause';
            else
                % Generate the ALGORITHM object
                try
                    item = obj.app.listD.items(1);
                    para = cell(1,length(item.edit));
                    for j = 1 : length(para)
                        if ~isempty(item.edit(j).Value)
                            para{j} = str2num(item.edit(j).Value);
                            assert(~isempty(para{j}),'The parameter "%s" of %s is illegal.',item.label(j).Text,item.title.Text);
                        end
                    end
                    ALG = feval(item.title.Text,'parameter',para(3:end),'outputFcn',@obj.outputFcn,'save',inf);
                catch err
                    uialert(obj.GUI.app.figure,err.message,'Invalid parameter settings');
                    return;
                end
                % Generate the PROBLEM object
                PRO = obj.cb_validation(para(1:2));
                % Update the data
                obj.data = {ALG,PRO};
                % Update the GUI
                [obj.GUI.app.button.Enable]   = deal('off');
                [obj.app.buttonA.Enable]      = deal('off');
                obj.app.editB(1).Enable       = 'off';
                [obj.app.editB(2:end).Enable] = deal('off');
                [obj.app.buttonB.Enable]      = deal('off');
                obj.app.dropB.Enable          = 'off';
                [obj.app.stateC.Enable]       = deal('off');
                obj.app.listC.Enable          = 'off';
                obj.app.listD.Enable          = 'off';
                [obj.app.toolD(1:2).Visible]  = deal('off');
                obj.app.buttonD(1).Text       = 'Pause';
                obj.app.buttonD(2).Enable     = 'on';
                obj.app.menu((PRO.M<=1)+1).value = 1;
                [obj.app.menu(1).items(3:end).Enable] = deal('off');
                [obj.app.menu(2).items(2:end).Enable] = deal('off');
                set([obj.pro.edit,obj.con.edit],'Enable','off');
                set([obj.pro.button,obj.con.button],'Enable','off');
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
            [obj.GUI.app.button.Enable]   = deal('on');
            [obj.app.buttonA.Enable]      = deal('on');
            obj.app.buttonA(2).Enable     = length(obj.pro) > 1;
            obj.app.buttonA(4).Enable     = ~isempty(obj.con);
            obj.app.editB(1).Enable       = 'on';
            [obj.app.editB(2:3).Enable]   = deal(obj.app.dropB.Value==1);
            [obj.app.editB(4:end).Enable] = deal('on');
            [obj.app.buttonB.Enable]      = deal('on');
            obj.app.dropB.Enable          = 'on';
            [obj.app.stateC.Enable]       = deal('on');
            obj.app.listC.Enable          = 'on';
            obj.app.listD.Enable          = 'on';
            [obj.app.toolD(1:2).Visible]  = deal('on');
            obj.app.buttonD(1).Text       = 'Start';
            obj.app.buttonD(2).Enable     = 'off';
            [obj.app.menu(1).items(3:end).Enable] = deal('on');
            [obj.app.menu(2).items(2:end).Enable] = deal('on');
            set([obj.pro.edit,obj.con.edit],'Enable','on');
            set([obj.pro.button,obj.con.button],'Enable','on');
            % Save the solutions to workspace
            assignin('base','X',obj.data{1}.result{end}.decs);
        end
        %% Output function
        function outputFcn(obj,Algorithm,Problem)
            obj.app.slider.Value = Problem.FE/max(Problem.FE,Problem.maxFE);
            obj.cb_slider();
            assert(strcmp(obj.app.buttonD(2).Enable,'on'),'PlatEMO:Termination','');
            if strcmp(obj.app.buttonD(1).Text,'Continue')
                waitfor(obj.app.buttonD(1),'Text');
            end
            assert(strcmp(obj.app.buttonD(2).Enable,'on'),'PlatEMO:Termination','');
        end
        %% Show the specified data
        function cb_slider(obj,~,~,ax)
            if ~isempty(obj.data)
                % Determine the current number of evaluationsnumber of evaluations
                ALG  = obj.data{1};
                PRO  = obj.data{2};
                rate = PRO.FE/max(PRO.FE,PRO.maxFE);
                obj.app.slider.Value      = min(obj.app.slider.Value,rate);
                obj.app.slider.MajorTicks = 0 : 0.25 : rate;
                obj.app.slider.MinorTicks = 0 : 0.01 : rate;
                index = max(1,round(obj.app.slider.Value/rate*size(ALG.result,1)));
                obj.app.labelD.Text = sprintf('%d evaluations',ALG.result{index,1});
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
                        otherwise
                            obj.app.waittip.Visible = 'on'; drawnow();
                            Draw(ALG.Metric(obj.app.menu(1).string,10),'-k.','LineWidth',1.5,'MarkerSize',10,{'Number of function evaluations',obj.app.menu(1).string,[]});
                            obj.app.waittip.Visible = 'off';
                    end
                else
                    % Show the result of single-objective optimization
                    switch obj.app.menu(2).value
                        case 1
                            PRO.DrawDec(ALG.result{index,2});
                        otherwise
                            Draw(ALG.Metric(obj.app.menu(2).string,10),'-k.','LineWidth',1.5,'MarkerSize',10,{'Number of function evaluations',obj.app.menu(2).string,[]});
                    end
                end
            end
        end
        %% Create the gif
        function cb_toolbutton1(obj,~,~)
            if ~isempty(obj.data)
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
            if ~isempty(obj.data)
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
            if ~isempty(obj.data)
                obj.app.menu((obj.data{2}.M<=1)+1).show();
            end
        end
    end
end

%% Get the function handle
function f = GetFcn(str,~)
    if ~isempty(regexp(str,'^<.+>$','once'))
        [folder,file] = fileparts(str(2:end-1));
        addpath(folder);
        f = str2func(file);
    elseif nargin == 1
        f = str2func(['@(x,data)',str]);
    else
        f = str2func(['@(n,data)',str]);
    end
end