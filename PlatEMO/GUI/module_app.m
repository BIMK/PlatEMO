classdef module_app < handle
%module_test - Application module.

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
            obj.app.grid(1)    = GUI.APP(2,1,uigridlayout(obj.app.maingrid,'RowHeight',{1,'1x',1},'ColumnWidth',{18,18,70,'1x'},'Padding',[5 5 5 5],'RowSpacing',0,'ColumnSpacing',7,'BackgroundColor',[.95 .95 1]));
            tempPanel          = GUI.APP(2,1,uipanel(obj.app.grid(1),'BorderType','none','BackgroundColor',[.95 .95 1]));
            obj.app.buttonA(1) = uibutton(tempPanel,'Position',[-2.5 -2.5 24 24],'Text','','Icon',obj.GUI.icon.loadtable,'BackgroundColor',[.95 .95 1],'Tooltip','Load a problem','Interruptible','off','BusyAction','cancel','ButtonpushedFcn',@obj.cb_loadProblem);
            tempPanel          = GUI.APP(2,2,uipanel(obj.app.grid(1),'BorderType','none','BackgroundColor',[.95 .95 1]));
            obj.app.buttonA(2) = uibutton(tempPanel,'Position',[-2.5 -2.5 24 24],'Text','','Icon',obj.GUI.icon.savetable,'BackgroundColor',[.95 .95 1],'Tooltip','Save the problem','Interruptible','off','BusyAction','cancel','ButtonpushedFcn',@obj.cb_saveProblem);
            obj.app.buttonA(3) = GUI.APP([1 3],3,uibutton(obj.app.grid(1),'Text','Validation','BackgroundColor',[.95 .95 1],'Tooltip','Check the validity of the problem','Interruptible','off','BusyAction','cancel','ButtonpushedFcn',@(~,~)obj.cb_validation()));
            obj.app.dropA      = GUI.APP([1 3],4,uidropdown(obj.app.grid(1),'BackgroundColor',[.95 .95 1],'Tooltip','Select predefined examples','Items',{'User-defined problem','Single-objective problem','Multi-objective problem','Many-objective problem','Constrained single-objective problem','Constrained multi-objective problem','Rotated optimization problem','Binary optimization problem','Permutation optimization problem'},'ItemsData',1:9,'Value',1,'Interruptible','off','BusyAction','cancel','ValueChangedFcn',@obj.cb_selectProblem));
          
            % The second panel
            obj.app.grid(2)     = GUI.APP(3,1,uigridlayout(obj.app.maingrid,'RowHeight',num2cell([25,25,25,25,25,25,25,0,25,0,25,0,25,25,zeros(1,19),25,25,zeros(1,19),25,25,25,25,25,25]),'ColumnWidth',{75,50,'1x',25,20},'Padding',[0 10 0 0],'RowSpacing',5,'ColumnSpacing',5,'Scrollable','on','BackgroundColor','w'));
            tempGrid            = GUI.APP(1,[1 5],uigridlayout(obj.app.grid(2),'RowHeight',{'1x'},'ColumnWidth',{'1x'},'Padding',[0 0 0 6],'RowSpacing',0,'ColumnSpacing',5,'BackgroundColor','w'));
            obj.app.titleB(1)   = GUI.APP(1,1,uilabel(tempGrid,'Text',' Decision vector','FontSize',13,'FontColor',[.9 .5 .2],'BackgroundColor',[.9 .9 .9],'FontWeight','bold'));
            obj.app.labelB(1)   = GUI.APP(2,1,uilabel(obj.app.grid(2),'Text','x =','HorizontalAlignment','right'));
            obj.app.editB(1)    = GUI.APP(2,2,uieditfield(obj.app.grid(2),'numeric','Value',10,'limits',[1 inf],'RoundFractionalValues','on'));
            obj.app.dropB       = GUI.APP(2,3,uidropdown(obj.app.grid(2),'BackgroundColor','w','Items',{'real','binary','permutation'},'ItemsData',1:3,'Value',1,'ValueChangedFcn',{@obj.cb_updateProblem,0}));
            obj.app.labelB(2)   = GUI.APP(2,[4 5],uilabel(obj.app.grid(2),'Text','variables'));
            obj.app.labelB(3)   = GUI.APP(3,1,uilabel(obj.app.grid(2),'Text','e.g.  x =','HorizontalAlignment','right'));
            obj.app.labelB(4)   = GUI.APP(3,[2 5],uilabel(obj.app.grid(2),'Text',['(',num2str(rand(1,10),'% .1f'),')']));
            tempGrid            = GUI.APP(4,[1 5],uigridlayout(obj.app.grid(2),'RowHeight',{'1x'},'ColumnWidth',{'1x'},'Padding',[0 0 0 6],'RowSpacing',0,'ColumnSpacing',5,'BackgroundColor','w'));
            obj.app.titleB(2)   = GUI.APP(1,1,uilabel(tempGrid,'Text',' Decision space','FontSize',13,'FontColor',[.9 .5 .2],'BackgroundColor',[.9 .9 .9],'FontWeight','bold'));
            obj.app.labelB(5)   = GUI.APP(5,1,uilabel(obj.app.grid(2),'Text','x >=','HorizontalAlignment','right'));
            obj.app.editB(2)    = GUI.APP(5,[2 4],uieditfield(obj.app.grid(2),'Value','0,0,0,0,0,0,0,0,0,0','Tooltip','Lower bound of each decision variable'));
            obj.app.labelB(6)   = GUI.APP(6,1,uilabel(obj.app.grid(2),'Text','x <=','HorizontalAlignment','right'));
            obj.app.editB(3)    = GUI.APP(6,[2 4],uieditfield(obj.app.grid(2),'Value','1,1,1,1,1,1,1,1,1,1','Tooltip','Upper bound of each decision variable'));
            tempGrid            = GUI.APP(7,[1 5],uigridlayout(obj.app.grid(2),'RowHeight',{'1x'},'ColumnWidth',{50,30,30,'1x'},'Padding',[0 0 0 6],'RowSpacing',0,'ColumnSpacing',5,'BackgroundColor','w'));
            obj.app.titleB(3)   = GUI.APP(1,[1 4],uilabel(tempGrid,'Text',' Dataset','FontSize',13,'FontColor',[.9 .5 .2],'BackgroundColor',[.9 .9 .9],'FontWeight','bold'));
            obj.app.buttonB(1)  = GUI.APP(1,2,uibutton(tempGrid,'Text','Add','FontSize',10,'Tooltip','Add a dataset','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_updateProblem,1}));
            obj.app.buttonB(2)  = GUI.APP(1,3,uibutton(tempGrid,'Text','Del','Enable','off','FontSize',10,'Tooltip','Delete the dataset','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_updateProblem,-1}));
            obj.app.labelB(7)   = GUI.APP(8,1,uilabel(obj.app.grid(2),'Text','data =','HorizontalAlignment','right'));
            obj.app.editB(4)    = GUI.APP(8,[2 4],uieditfield(obj.app.grid(2),'Value','eye(10)','Tooltip','User-defined parameters of the problem'));
            obj.app.buttonB(3)  = GUI.APP(8,5,uibutton(obj.app.grid(2),'Text','...','Tooltip','Load from file','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_loadFunction,obj.app.editB(4),{'*.mat';'*.txt'}}));
            tempGrid            = GUI.APP(9,[1 5],uigridlayout(obj.app.grid(2),'RowHeight',{'1x'},'ColumnWidth',{133,30,30,'1x'},'Padding',[0 0 0 6],'RowSpacing',0,'ColumnSpacing',5,'BackgroundColor','w'));
            obj.app.titleB(4)   = GUI.APP(1,[1 4],uilabel(tempGrid,'Text',' Initialization function','FontSize',13,'FontColor',[.9 .5 .2],'BackgroundColor',[.9 .9 .9],'FontWeight','bold'));
            obj.app.buttonB(4)  = GUI.APP(1,2,uibutton(tempGrid,'Text','Add','FontSize',10,'Tooltip','Add an initialization function','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_updateProblem,2}));
            obj.app.buttonB(5)  = GUI.APP(1,3,uibutton(tempGrid,'Text','Del','Enable','off','FontSize',10,'Tooltip','Delete the initialization function','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_updateProblem,-2}));
            obj.app.labelB(8)   = GUI.APP(10,1,uilabel(obj.app.grid(2),'Text','I(N) =','HorizontalAlignment','right'));
            obj.app.editB(5)    = GUI.APP(10,[2 4],uieditfield(obj.app.grid(2),'Value','rand(N,10)','Tooltip','Function for generating an initial population with size n'));
            obj.app.buttonB(6)  = GUI.APP(10,5,uibutton(obj.app.grid(2),'Text','...','Tooltip','Load from file','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_loadFunction,obj.app.editB(5),'*.m'}));
            tempGrid            = GUI.APP(11,[1 5],uigridlayout(obj.app.grid(2),'RowHeight',{'1x'},'ColumnWidth',{100,30,30,'1x'},'Padding',[0 0 0 6],'RowSpacing',0,'ColumnSpacing',5,'BackgroundColor','w'));
            obj.app.titleB(5)   = GUI.APP(1,[1 4],uilabel(tempGrid,'Text',' Repair function','FontSize',13,'FontColor',[.9 .5 .2],'BackgroundColor',[.9 .9 .9],'FontWeight','bold'));
            obj.app.buttonB(7)  = GUI.APP(1,2,uibutton(tempGrid,'Text','Add','FontSize',10,'Tooltip','Add a repair function','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_updateProblem,3}));
            obj.app.buttonB(8)  = GUI.APP(1,3,uibutton(tempGrid,'Text','Del','Enable','off','FontSize',10,'Tooltip','Delete the repair function','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_updateProblem,-3}));
            obj.app.labelB(9)   = GUI.APP(12,1,uilabel(obj.app.grid(2),'Text','R(x) =','HorizontalAlignment','right'));
            obj.app.editB(6)    = GUI.APP(12,[2 4],uieditfield(obj.app.grid(2),'Value','max(min(x,1),0)','Tooltip','Function for repairing invalid solution'));
            obj.app.buttonB(9)  = GUI.APP(12,5,uibutton(obj.app.grid(2),'Text','...','Tooltip','Load from file','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_loadFunction,obj.app.editB(6),'*.m'}));
            tempGrid            = GUI.APP(13,[1 5],uigridlayout(obj.app.grid(2),'RowHeight',{'1x'},'ColumnWidth',{123,30,30,'1x'},'Padding',[0 0 0 6],'RowSpacing',0,'ColumnSpacing',5,'BackgroundColor','w'));
            obj.app.titleB(6)   = GUI.APP(1,[1 4],uilabel(tempGrid,'Text',' Objective functions','FontSize',13,'FontColor',[.9 .5 .2],'BackgroundColor',[.9 .9 .9],'FontWeight','bold'));
            obj.app.buttonB(10) = GUI.APP(1,2,uibutton(tempGrid,'Text','Add','FontSize',10,'Tooltip','Add an objective function','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_updateProblem,4}));
            obj.app.buttonB(11) = GUI.APP(1,3,uibutton(tempGrid,'Text','Del','Enable','off','FontSize',10,'Tooltip','Delete the last objective function','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_updateProblem,-4}));
            obj.pro(1).label    = GUI.APP(14,1,uilabel(obj.app.grid(2),'Text','f1(x) =','HorizontalAlignment','right'));
            obj.pro(1).edit     = GUI.APP(14,[2 4],uieditfield(obj.app.grid(2),'Value','mean(x)','Tooltip','Objective function to be minimized'));
            obj.pro(1).button   = GUI.APP(14,5,uibutton(obj.app.grid(2),'Text','...','Tooltip','Load from file','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_loadFunction,obj.pro(1).edit,'*.m'}));
            tempGrid            = GUI.APP(34,[1 5],uigridlayout(obj.app.grid(2),'RowHeight',{'1x'},'ColumnWidth',{130,30,30,'1x'},'Padding',[0 0 0 6],'RowSpacing',0,'ColumnSpacing',5,'BackgroundColor','w'));
            obj.app.titleB(7)   = GUI.APP(1,[1 4],uilabel(tempGrid,'Text',' Constraint functions','FontSize',13,'FontColor',[.9 .5 .2],'BackgroundColor',[.9 .9 .9],'FontWeight','bold'));
            obj.app.buttonB(12) = GUI.APP(1,2,uibutton(tempGrid,'Text','Add','FontSize',10,'Tooltip','Add a constraint function','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_updateProblem,5}));
            obj.app.buttonB(13) = GUI.APP(1,3,uibutton(tempGrid,'Text','Del','FontSize',10,'Tooltip','Delete the last constraint function','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_updateProblem,-5}));
            obj.con(1).labelA   = GUI.APP(35,1,uilabel(obj.app.grid(2),'Text','g1(x) =','HorizontalAlignment','right'));
            obj.con(1).edit     = GUI.APP(35,[2 3],uieditfield(obj.app.grid(2),'Value','0.5-mean(x)','Tooltip','Constraint function to be satisfied'));
            obj.con(1).labelB   = GUI.APP(35,4,uilabel(obj.app.grid(2),'Text','<= 0'));
            obj.con(1).button   = GUI.APP(35,5,uibutton(obj.app.grid(2),'Text','...','Tooltip','Load from file','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_loadFunction,obj.con(1).edit,'*.m'}));
            tempGrid            = GUI.APP(55,[1 5],uigridlayout(obj.app.grid(2),'RowHeight',{'1x'},'ColumnWidth',{127,30,30,'1x'},'Padding',[0 0 0 6],'RowSpacing',0,'ColumnSpacing',5,'BackgroundColor','w'));
            obj.app.titleB(7)   = GUI.APP(1,[1 4],uilabel(tempGrid,'Text',' Special difficulties','FontSize',13,'FontColor',[.9 .5 .2],'BackgroundColor',[.9 .9 .9],'FontWeight','bold'));
            obj.app.checkB(1)   = GUI.APP(56,[1 5],uicheckbox(obj.app.grid(2),'FontSize',11,'Text','<large> The search space is large, e.g., having more than 100 real decision variables.','WordWrap','on','ValueChangedFcn',@obj.cb_updateFilter));
            obj.app.checkB(2)   = GUI.APP(57,[1 5],uicheckbox(obj.app.grid(2),'FontSize',11,'Text','<expensive> The functions are computationally expensive, e.g., using less than 1000 function evaluations.','WordWrap','on','ValueChangedFcn',@obj.cb_updateFilter));
            obj.app.checkB(3)   = GUI.APP(58,[1 5],uicheckbox(obj.app.grid(2),'FontSize',11,'Text','<multimodal> There are multiple optimal solutions with similar objective values but different decision vectors, all of which should be obtained.','WordWrap','on','ValueChangedFcn',@obj.cb_updateFilter));
            obj.app.checkB(4)   = GUI.APP(59,[1 5],uicheckbox(obj.app.grid(2),'FontSize',11,'Text','<sparse> The optimal solutions are sparse, i.e., most decision variables of which are zero.','WordWrap','on','ValueChangedFcn',@obj.cb_updateFilter));
            obj.app.checkB(5)   = GUI.APP(60,[1 5],uicheckbox(obj.app.grid(2),'FontSize',11,'Text','<preference> Only searching for preferred regions of the Pareto front.','WordWrap','on','ValueChangedFcn',@obj.cb_updateFilter));
            
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
            obj.app.listD      = uilist(obj.app.grid(3),obj.GUI.app.figure,obj.GUI.icon);
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
            obj.app.toolD(1)   = axtoolbarbtn(tempTb,'push','Icon',obj.GUI.icon.gif,'Tooltip','Save the evolutionary process to gif','ButtonPushedFcn',@obj.cb_toolbutton1);
            obj.app.toolD(2)   = axtoolbarbtn(tempTb,'push','Icon',obj.GUI.icon.newfigure,'Tooltip','Open in new figure and save to workspace','ButtonPushedFcn',@obj.cb_toolbutton2);
            obj.app.toolD(3)   = axtoolbarbtn(tempTb,'push','Icon',obj.GUI.icon.datasource,'Tooltip','Data source','ButtonPushedFcn',@obj.cb_toolbutton3);
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
        %% Update the problem definition panel
        function cb_updateProblem(obj,~,~,index,noFlush)
            switch index
                case 0          % Encoding
                    obj.app.grid(2).RowHeight{4} = 25*(obj.app.dropB.Value<=1);
                	obj.app.grid(2).RowHeight{5} = 25*(obj.app.dropB.Value<=1);
                    obj.app.grid(2).RowHeight{6} = 25*(obj.app.dropB.Value<=1);
                    switch obj.app.dropB.Value
                        case 1
                            obj.app.labelB(4).Text = ['(',num2str(rand(1,obj.app.editB(1).Value),'% .1f'),')'];
                        case 2
                            obj.app.labelB(4).Text = ['(',num2str(rand(1,obj.app.editB(1).Value)<0.5,'% d'),')'];
                        case 3
                            obj.app.labelB(4).Text = ['(',num2str(randperm(obj.app.editB(1).Value),'% d'),')'];
                    end
                case {1,-1}     % Dataset
                    obj.app.grid(2).RowHeight{8} = max(0,25*sign(index));
                    obj.app.buttonB(1).Enable    = index < 0;
                    obj.app.buttonB(2).Enable    = index > 0;
                    if index > 0
                        str = ',data';
                    else
                        str = '';
                    end
                    obj.app.labelB(8).Text = sprintf('I(N%s) =',str);
                    obj.app.labelB(9).Text = sprintf('R(x%s) =',str);
                    for i = 1 : length(obj.pro)
                        obj.pro(i).label.Text = sprintf('f%d(x%s) =',i,str);
                    end
                    for i = 1 : length(obj.con)
                        obj.con(i).labelA.Text = sprintf('g%d(x%s) =',i,str);
                    end
                case {2,-2}     % Initialization function
                    obj.app.grid(2).RowHeight{10} = max(0,25*sign(index));
                    obj.app.buttonB(4).Enable     = index < 0;
                    obj.app.buttonB(5).Enable     = index > 0;
                case {3,-3}     % Repair function
                    obj.app.grid(2).RowHeight{12} = max(0,25*sign(index));
                    obj.app.buttonB(7).Enable     = index < 0;
                    obj.app.buttonB(8).Enable     = index > 0;
                case {4,-4}     % Objective functions
                    if index>0 && length(obj.pro)<20
                        if obj.app.grid(2).RowHeight{8} > 0
                            str = {',data','*data'};
                        else
                            str = {'',''};
                        end
                        item.label  = GUI.APP(length(obj.pro)+14,1,uilabel(obj.app.grid(2),'Text',sprintf('f%d(x%s) =',length(obj.pro)+1,str{1}),'HorizontalAlignment','right'));
                        item.edit   = GUI.APP(length(obj.pro)+14,[2 4],uieditfield(obj.app.grid(2),'Value',sprintf('mean(x%s)',str{2}),'Tooltip','Objective function to be minimized'));
                        item.button = GUI.APP(length(obj.pro)+14,5,uibutton(obj.app.grid(2),'Text','...','Tooltip','Load from file','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_loadFunction,item.edit,'*.m'}));
                        obj.pro     = [obj.pro,item];
                    elseif index<0 && length(obj.pro)>0
                        delete(obj.pro(end).label);
                        delete(obj.pro(end).edit);
                        delete(obj.pro(end).button);
                        obj.pro(end) = [];
                    end
                    obj.app.grid(2).RowHeight(14:33) = num2cell([zeros(1,length(obj.pro))+25,zeros(1,20-length(obj.pro))]);
                    obj.app.buttonB(10).Enable       = length(obj.pro) < 20;
                    obj.app.buttonB(11).Enable       = length(obj.pro) > 1;
                case {5,-5}     % Constraint functions
                    if index>0 && length(obj.con)<20
                        if obj.app.grid(2).RowHeight{8} > 0
                            str = {',data','*data'};
                        else
                            str = {'',''};
                        end
                        item.labelA = GUI.APP(length(obj.con)+35,1,uilabel(obj.app.grid(2),'Text',sprintf('g%d(x%s) =',length(obj.con)+1,str{1}),'HorizontalAlignment','right'));
                        item.edit   = GUI.APP(length(obj.con)+35,[2 3],uieditfield(obj.app.grid(2),'Value',sprintf('0.5-mean(x%s)',str{2}),'Tooltip','Constraint function to be satisfied'));
                        item.labelB = GUI.APP(length(obj.con)+35,4,uilabel(obj.app.grid(2),'Text','<= 0'));
                        item.button = GUI.APP(length(obj.con)+35,5,uibutton(obj.app.grid(2),'Text','...','Tooltip','Load from file','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_loadFunction,item.edit,'*.m'}));
                        obj.con     = [obj.con,item];
                    elseif index<0 && length(obj.con)>0
                        delete(obj.con(end).labelA);
                        delete(obj.con(end).edit);
                        delete(obj.con(end).labelB);
                        delete(obj.con(end).button);
                        obj.con(end) = [];
                    end
                    obj.app.grid(2).RowHeight(35:54) = num2cell([zeros(1,length(obj.con))+25,zeros(1,20-length(obj.con))]);
                    obj.app.buttonB(12).Enable       = length(obj.con) < 20;
                    obj.app.buttonB(13).Enable       = length(obj.con) > 0;
            end
            if nargin < 5 || ~noFlush
                obj.cb_updateFilter();
            end
        end
        %% Validate the problem
        function PRO = cb_validation(obj,para)
            if nargin < 2; para = {1,1}; end
            Str = GetSetProblem(obj);
            try
                PRO = UserProblem('N',para{1},'D',Str{1},'maxFE',para{2},'encoding',Str{2},'lower',Str{3},'upper',Str{4},'parameter',Str{5},'initFcn',Str{6},'decFcn',Str{7},'objFcn',Str{8},'conFcn',Str{9});
                PROBLEM.Current(PRO);
                PRO.Initialization();
            catch err
                if isempty(err.cause)
                    err = addCause(err,MException('','The problem definition is invalid'));
                end
                uialert(obj.GUI.app.figure,sprintf('%s, since %s',err.cause{1}.message,err.message),'Invalid problem');
                rethrow(err);
            end
            if nargout < 1
                uialert(obj.GUI.app.figure,'The problem definition is valid.','Valid problem','Icon','success');
            end
        end
        %% Define a problem
        function cb_selectProblem(obj,~,~)
            if obj.app.dropA.Value == 1
                if ~isempty(obj.app.dropA.UserData)
                    GetSetProblem(obj,obj.app.dropA.UserData);
                    obj.app.dropA.UserData = [];
                end
            else
                switch obj.app.dropA.Value
                    case 2
                        Str = {10,'real','zeros(1,10)','ones(1,10)','','','','mean(x.^2)',''};
                    case 3
                        Str = {15,'real','zeros(1,15)','ones(1,15)','','','',{'x(1)+mean(x(2:end).^2)','1-x(1).^2+mean(x(2:end).^2)'},''};
                    case 4
                        Str = {20,'real','zeros(1,20)','ones(1,20)','','','',{'x(1)*x(2)*x(3)+mean(x(4:end).^2)','x(1)*x(2)*(1-x(3))+mean(x(4:end).^2)','x(1)*(1-x(2))+mean(x(4:end).^2)','1-x(1)+mean(x(4:end).^2)'},''};
                    case 5
                        Str = {25,'real','zeros(1,25)','ones(1,25)','','','','mean(x.^2)','0.2-mean(x.^2)'};
                    case 6
                        Str = {30,'real','zeros(1,30)','ones(1,30)','','','',{'x(1)+mean(x(2:end).^2)','1-x(1).^2+mean(x(2:end).^2)'},{'0.2-mean(x(2:end).^2)','mean(x(2:end).^2)-0.8'}};
                    case 7
                        Str = {35,'real','zeros(1,35)','ones(1,35)','rand(35)','','','mean((x*data).^2)',''};
                    case 8
                        Str = {100,'binary','','','','rand(N,100)<0.2','','1-mean(x)','mean(x)-0.1'};
                    case 9
                        Str = {20,'permutation','','','','','','mean(cumsum(movsum(x,3)))',''};
                end
                if isempty(obj.app.dropA.UserData)
                    obj.app.dropA.UserData = GetSetProblem(obj);
                end
                GetSetProblem(obj,Str);
            end
        end
        %% Load a problem
        function cb_loadProblem(obj,~,~)
            [file,folder] = uigetfile('*.m');
            if ischar(file)
                try
                    f = fopen(fullfile(folder,file));
                    fgetl(f); fgetl(f); fgetl(f); fgetl(f); fgetl(f); fgetl(f);
                    Str{1} = eval(regexp(fgetl(f),'(?<=''D'',).*(?=,...)','match','once'));
                    Str{2} = regexp(fgetl(f),'(?<=''encoding'','').*(?='',...)','match','once');
                    Str{3} = regexp(fgetl(f),'(?<=''lower'','').*(?='',...)','match','once');
                    Str{4} = regexp(fgetl(f),'(?<=''upper'','').*(?='',...)','match','once');
                    Str{5} = regexp(fgetl(f),'(?<=''parameter'','').*(?='',...)','match','once');
                    Str{6} = regexp(fgetl(f),'(?<=''initFcn'','').*(?='',...)','match','once');
                    Str{7} = regexp(fgetl(f),'(?<=''decFcn'','').*(?='',...)','match','once');
                    Str{8} = eval(regexp(fgetl(f),'(?<=''objFcn'',).*(?=,...)','match','once'));
                    Str{9} = eval(regexp(fgetl(f),'(?<=''conFcn'',).*(?=\);)','match','once'));
                    obj.app.dropA.Value = 1;
                    GetSetProblem(obj,Str);
                catch err
                    uialert(obj.GUI.app.figure,'Fail to load the problem, please refer to the command window for details.','Error');
                    rethrow(err);
                end
            end
        end
        %% Save the problem
        function cb_saveProblem(obj,~,~)
            [Name,Path] = uiputfile('*.m','','new');
            if ischar(Name)
                try
                    Str      = GetSetProblem(obj);
                    [~,name] = fileparts(Name);
                    Code = {['classdef ',name,' < UserProblem % < PROBLEM']
                           '% The code is automatically generated and cannot be modified'
                           ''
                           '    methods'
                           ['        function obj = ',name,'(varargin)']
                           ['            obj = obj@UserProblem(varargin{:},...']
                           ['            ''D'',',num2str(Str{1}),',...']
                           ['            ''encoding'',''',Str{2},''',...']
                           ['            ''lower'',''',Str{3},''',...']
                           ['            ''upper'',''',Str{4},''',...']
                           ['            ''parameter'',''',Str{5},''',...']
                           ['            ''initFcn'',''',Str{6},''',...']
                           ['            ''decFcn'',''',Str{7},''',...']
                           ['            ''objFcn'',{',strjoin(strcat('''',Str{8},''''),','),'},...']
                           ['            ''conFcn'',{',strjoin(strcat('''',Str{9},''''),','),'});']
                           '        end'
                           '    end'
                           'end'};
                    fid = fopen(fullfile(Path,Name),'wt');
                    for i = 1 : length(Code)
                        fprintf(fid,'%s\n',Code{i});
                    end
                    fclose(fid);
                catch err
                    uialert(obj.GUI.app.figure,'Fail to save the problem, please refer to the command window for details.','Error');
                    rethrow(err);
                end
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
            obj.app.stateC(7).Value  = obj.app.checkB(1).Value;
            obj.app.stateC(8).Value  = ~isempty(obj.con);
            obj.app.stateC(9).Value  = obj.app.checkB(2).Value;
            obj.app.stateC(10).Value = obj.app.checkB(3).Value;
            obj.app.stateC(11).Value = obj.app.checkB(4).Value;
            obj.app.stateC(12).Value = obj.app.checkB(5).Value;
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
                obj.data = {ALG,PRO,{obj.app.buttonB.Enable}};
                % Update the GUI
                [obj.GUI.app.button.Enable]   = deal('off');
                [obj.app.buttonA.Enable]      = deal('off');
                obj.app.dropA.Enable          = 'off';
                obj.app.editB(1).Enable       = 'off';
                [obj.app.editB(2:end).Enable] = deal('off');
                [obj.app.buttonB.Enable]      = deal('off');
                [obj.app.checkB.Enable]       = deal('off');
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
            obj.app.dropA.Enable          = 'on';
            obj.app.editB(1).Enable       = 'on';
            [obj.app.editB(2:end).Enable] = deal('on');
            [obj.app.buttonB.Enable]      = deal(obj.data{3}{:});
            [obj.app.checkB.Enable]       = deal('on');
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

function Str = GetSetProblem(obj,Str)
    if nargout > 0  % Get the problem
        Str{1} = obj.app.editB(1).Value;                    % No. of variable
        Str{2} = obj.app.dropB.Items{obj.app.dropB.Value};  % Encoding
        Str{3} = obj.app.editB(2).Value;                    % Lower bound
        Str{4} = obj.app.editB(3).Value;                    % Upper bound
        if obj.app.grid(2).RowHeight{8} > 0                 % Dataset
            Str{5} = obj.app.editB(4).Value;
        else
            Str{5} = '';
        end
        if obj.app.grid(2).RowHeight{10} > 0                % Initialization function
            Str{6} = obj.app.editB(5).Value;
        else
            Str{6} = '';
        end
        if obj.app.grid(2).RowHeight{12} > 0                % Repair function
            Str{7} = obj.app.editB(6).Value;
        else
            Str{7} = '';
        end
        Str{8} = get([obj.pro.edit],'Value');               % Objective functions
        if ~iscell(Str{8}) && ~isempty(Str{8}); Str{8} = Str(8); end
        Str{9} = get([obj.con.edit],'Value');               % Constraint functions
        if ~iscell(Str{9}) && ~isempty(Str{9}); Str{9} = Str(9); end
    end
    if nargin > 1	% Set the problem
        obj.app.editB(1).Value  = Str{1};                               % Number of decision variables
        [~,obj.app.dropB.Value] = ismember(Str{2},obj.app.dropB.Items);	% Encoding
        obj.cb_updateProblem([],[],0,1);
        obj.app.editB(2).Value  = Str{3};                               % Lower bound
        obj.app.editB(3).Value  = Str{4};                             	% Upper bound
        obj.app.editB(4).Value  = Str{5};                              	% Dataset
        obj.cb_updateProblem([],[],1*(-1)^isempty(obj.app.editB(4).Value),1);
        obj.app.editB(5).Value  = Str{6};                            	% Initialization function
        obj.cb_updateProblem([],[],2*(-1)^isempty(obj.app.editB(5).Value),1);
        obj.app.editB(6).Value  = Str{7};                              	% Repair function
        obj.cb_updateProblem([],[],3*(-1)^isempty(obj.app.editB(6).Value),1);
        for i = 1 : length(obj.pro)                                     % Objective functions
            obj.cb_updateProblem([],[],-4,1);
        end
        if ~iscell(Str{8}); Str{8} = Str(8); end
        Str{8}(cellfun(@isempty,Str{8})) = [];
        for i = 1 : length(Str{8})
            obj.cb_updateProblem([],[],4,1);
            obj.pro(end).edit.Value = Str{8}{i};
        end
        for i = 1 : length(obj.con)                                     % Constraint functions
            obj.cb_updateProblem([],[],-5,1);
        end
        if ~iscell(Str{9}); Str{9} = Str(9); end
        Str{9}(cellfun(@isempty,Str{9})) = [];
        for i = 1 : length(Str{9})
            obj.cb_updateProblem([],[],5,1);
            obj.con(end).edit.Value = Str{9}{i};
        end
        obj.cb_updateFilter();
    end
end