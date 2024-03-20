classdef module_app < handle
%module_test - Application module.

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
        pro  = [];          % Components of problem functions
        con  = [];          % Components of constraint functions
        data = {};          % All the results
    end
    methods(Access = ?GUI)
        %% Constructor
        function obj = module_app(GUI)
            % The main grid
            obj.GUI = GUI;
            obj.app.maingrid = GUI.APP(3,1,uigridlayout(obj.GUI.app.maingrid,'RowHeight',{20,30,'1x'},'ColumnWidth',{'2x',5,1,'1.2x',1,'2x','0.8x','0.2x'},'Padding',[5 5 5 5],'RowSpacing',5,'ColumnSpacing',0,'BackgroundColor','w'));
            obj.app.label(1) = GUI.APP(1,1,uilabel(obj.app.maingrid,'Text','Problem definition','HorizontalAlignment','center','FontSize',11,'FontColor',[.5 .5 .5]));
            obj.app.label(2) = GUI.APP(1,4,uilabel(obj.app.maingrid,'Text','Algorithm selection','HorizontalAlignment','center','FontSize',11,'FontColor',[.5 .5 .5]));
            obj.app.label(3) = GUI.APP(1,[6 8],uilabel(obj.app.maingrid,'Text','Result display','HorizontalAlignment','center','FontSize',11,'FontColor',[.5 .5 .5]));
            GUI.APP([1 3],3,uipanel(obj.app.maingrid,'BackgroundColor',[.8 .8 .8]));
            GUI.APP([1 3],5,uipanel(obj.app.maingrid,'BackgroundColor',[.8 .8 .8]));
            
            % The first panel
            obj.app.grid(1)    = GUI.APP(2,1,uigridlayout(obj.app.maingrid,'RowHeight',{1,'1x',1},'ColumnWidth',{18,18,70,'1x'},'Padding',[5 5 5 5],'RowSpacing',0,'ColumnSpacing',7,'BackgroundColor',[.95 .95 1]));
            tempPanel          = GUI.APP(2,1,uipanel(obj.app.grid(1),'BorderType','none','BackgroundColor',[.95 .95 1]));
            obj.app.buttonA(1) = uibutton(tempPanel,'Position',[-2.5 -2.5 24 24],'Text','','Icon',obj.GUI.icon.loadtable,'BackgroundColor',[.95 .95 1],'Tooltip','Load a problem','Interruptible','off','BusyAction','cancel','ButtonpushedFcn',@obj.cb_loadProblem);
            tempPanel          = GUI.APP(2,2,uipanel(obj.app.grid(1),'BorderType','none','BackgroundColor',[.95 .95 1]));
            obj.app.buttonA(2) = uibutton(tempPanel,'Position',[-2.5 -2.5 24 24],'Text','','Icon',obj.GUI.icon.savetable,'BackgroundColor',[.95 .95 1],'Tooltip','Save the problem','Interruptible','off','BusyAction','cancel','ButtonpushedFcn',@obj.cb_saveProblem);
            obj.app.buttonA(3) = GUI.APP([1 3],3,uibutton(obj.app.grid(1),'Text','Validation','BackgroundColor',[.95 .95 1],'Tooltip','Check the validity of the problem','Interruptible','off','BusyAction','cancel','ButtonpushedFcn',@(~,~)obj.cb_validation()));
            obj.app.dropA      = GUI.APP([1 3],4,uidropdown(obj.app.grid(1),'BackgroundColor',[.95 .95 1],'Tooltip','Select predefined examples','Items',{'User-defined problem','Single-objective optimization','Multi-objective optimization','Many-objective optimization','Constrained single-objective optimization',...
                                         'Constrained multi-objective optimization','Rotated optimization','Integer optimization','Binary optimization','Permutation optimization','Hybrid optimization'},'ItemsData',1:11,'Value',1,'Interruptible','off','BusyAction','cancel','ValueChangedFcn',@obj.cb_selectProblem));
          
            % The second panel
            obj.app.grid(2)     = GUI.APP(3,1,uigridlayout(obj.app.maingrid,'RowHeight',num2cell([25,25,25,25,25,25,25,0,25,0,25,0,25,25,zeros(1,19),25,25,zeros(1,19),25,25,25,25,25,25]),'ColumnWidth',{75,50,'1x',25,20,20},'Padding',[0 10 0 0],'RowSpacing',5,'ColumnSpacing',5,'Scrollable','on','BackgroundColor','w'));
            % Encoding scheme
            tempGrid            = GUI.APP(1,[1 6],uigridlayout(obj.app.grid(2),'RowHeight',{'1x'},'ColumnWidth',{'1x'},'Padding',[0 0 0 6],'RowSpacing',0,'ColumnSpacing',5,'BackgroundColor','w'));
            obj.app.titleB(1)   = GUI.APP(1,1,uilabel(tempGrid,'Text',' Encoding scheme','FontSize',13,'FontColor',[.9 .5 .2],'BackgroundColor',[.9 .9 .9],'FontWeight','bold'));
            obj.app.labelB(1)   = GUI.APP(2,1,uilabel(obj.app.grid(2),'Text','Encoding: x =','HorizontalAlignment','right'));
            obj.app.editB(1)    = GUI.APP(2,[2 5],uieditfield(obj.app.grid(2),'Value','1,1,1,1,1,1,1,1,1,1','Tooltip','Type of each decision variable','ValueChangedFcn',@obj.cb_updateFilter));
            obj.app.labelB(3)   = GUI.APP(3,[1 6],uilabel(obj.app.grid(2),'Text','(1. real number 2. integer 3. label 4. binary number 5. permutation)'));
            % Decision space
            tempGrid            = GUI.APP(4,[1 6],uigridlayout(obj.app.grid(2),'RowHeight',{'1x'},'ColumnWidth',{'1x'},'Padding',[0 0 0 6],'RowSpacing',0,'ColumnSpacing',5,'BackgroundColor','w'));
            obj.app.titleB(2)   = GUI.APP(1,1,uilabel(tempGrid,'Text',' Decision space','FontSize',13,'FontColor',[.9 .5 .2],'BackgroundColor',[.9 .9 .9],'FontWeight','bold'));
            obj.app.labelB(4)   = GUI.APP(5,1,uilabel(obj.app.grid(2),'Text','Lower: x >=','HorizontalAlignment','right'));
            obj.app.editB(2)    = GUI.APP(5,[2 5],uieditfield(obj.app.grid(2),'Value','0,0,0,0,0,0,0,0,0,0','Tooltip','Lower bound of each decision variable'));
            obj.app.labelB(5)   = GUI.APP(6,1,uilabel(obj.app.grid(2),'Text','Upper: x <=','HorizontalAlignment','right'));
            obj.app.editB(3)    = GUI.APP(6,[2 5],uieditfield(obj.app.grid(2),'Value','1,1,1,1,1,1,1,1,1,1','Tooltip','Upper bound of each decision variable'));
            % Data
            tempGrid            = GUI.APP(7,[1 6],uigridlayout(obj.app.grid(2),'RowHeight',{'1x'},'ColumnWidth',{50,30,30,'1x'},'Padding',[0 0 0 6],'RowSpacing',0,'ColumnSpacing',5,'BackgroundColor','w'));
            obj.app.titleB(3)   = GUI.APP(1,[1 4],uilabel(tempGrid,'Text',' Data','FontSize',13,'FontColor',[.9 .5 .2],'BackgroundColor',[.9 .9 .9],'FontWeight','bold'));
            obj.app.buttonB(1)  = GUI.APP(1,2,uibutton(tempGrid,'Text','Add','FontSize',10,'Tooltip','Add a data','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_updateProblem,1}));
            obj.app.buttonB(2)  = GUI.APP(1,3,uibutton(tempGrid,'Text','Del','Enable',false,'FontSize',10,'Tooltip','Delete the data','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_updateProblem,-1}));
            obj.app.labelB(6)   = GUI.APP(8,1,uilabel(obj.app.grid(2),'Text','data =','HorizontalAlignment','right'));
            obj.app.editB(4)    = GUI.APP(8,[2 4],uieditfield(obj.app.grid(2),'Value','eye(10)','Tooltip','User-defined data of the problem'));
            obj.app.buttonB(3)  = GUI.APP(8,5,uibutton(obj.app.grid(2),'Text','...','Tooltip','Load existing data','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_loadFunction,obj.app.editB(4),{'*.txt;*.dat;*.csv','Text file';'*.mat','MAT file'}}));
            % Initialization function
            tempGrid            = GUI.APP(9,[1 6],uigridlayout(obj.app.grid(2),'RowHeight',{'1x'},'ColumnWidth',{135,30,30,'1x'},'Padding',[0 0 0 6],'RowSpacing',0,'ColumnSpacing',5,'BackgroundColor','w'));
            obj.app.titleB(4)   = GUI.APP(1,[1 4],uilabel(tempGrid,'Text',' Initialization function','FontSize',13,'FontColor',[.9 .5 .2],'BackgroundColor',[.9 .9 .9],'FontWeight','bold'));
            obj.app.buttonB(4)  = GUI.APP(1,2,uibutton(tempGrid,'Text','Add','FontSize',10,'Tooltip','Add an initialization function','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_updateProblem,2}));
            obj.app.buttonB(5)  = GUI.APP(1,3,uibutton(tempGrid,'Text','Del','Enable',false,'FontSize',10,'Tooltip','Delete the initialization function','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_updateProblem,-2}));
            obj.app.labelB(7)   = GUI.APP(10,1,uilabel(obj.app.grid(2),'Text','I(N) =','HorizontalAlignment','right'));
            obj.app.editB(5)    = GUI.APP(10,[2 4],uieditfield(obj.app.grid(2),'Value','rand(N,10)','Tooltip','Function for generating an initial population with size N'));
            obj.app.buttonB(6)  = GUI.APP(10,5,uibutton(obj.app.grid(2),'Text','...','Tooltip','Load existing function','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_loadFunction,obj.app.editB(5),{'*.m','MATLAB function'}}));
            obj.app.buttonB(7)  = GUI.APP(10,6,uibutton(obj.app.grid(2),'Text','+','Fontsize',15,'Tooltip','Create new function','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_createFunction,obj.app.editB(5),'initialization'}));
            % Repair function
            tempGrid            = GUI.APP(11,[1 6],uigridlayout(obj.app.grid(2),'RowHeight',{'1x'},'ColumnWidth',{135,30,30,50,'1x'},'Padding',[0 0 0 6],'RowSpacing',0,'ColumnSpacing',5,'BackgroundColor','w'));
            obj.app.titleB(5)   = GUI.APP(1,[1 5],uilabel(tempGrid,'Text',' Repair function','FontSize',13,'FontColor',[.9 .5 .2],'BackgroundColor',[.9 .9 .9],'FontWeight','bold'));
            obj.app.buttonB(8)  = GUI.APP(1,2,uibutton(tempGrid,'Text','Add','FontSize',10,'Tooltip','Add a repair function','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_updateProblem,3}));
            obj.app.buttonB(9)  = GUI.APP(1,3,uibutton(tempGrid,'Text','Del','Enable',false,'FontSize',10,'Tooltip','Delete the repair function','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_updateProblem,-3}));
            obj.app.buttonB(10) = GUI.APP(1,4,uibutton(tempGrid,'Text','Integrate','FontSize',10,'Tooltip','Integrate the repair function, objective functions, and constraint functions into a single evaluation function','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_updateProblem,6}));
            obj.app.labelB(8)   = GUI.APP(12,1,uilabel(obj.app.grid(2),'Text','R(x) =','HorizontalAlignment','right'));
            obj.app.editB(6)    = GUI.APP(12,[2 4],uieditfield(obj.app.grid(2),'Value','max(min(x,1),0)','Tooltip','Function for repairing invalid solution','UserData',''));
            obj.app.buttonB(11) = GUI.APP(12,5,uibutton(obj.app.grid(2),'Text','...','Tooltip','Load existing function','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_loadFunction,obj.app.editB(6),{'*.m','MATLAB function'}}));
            obj.app.buttonB(12) = GUI.APP(12,6,uibutton(obj.app.grid(2),'Text','+','Fontsize',15,'Tooltip','Create new function','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_createFunction,obj.app.editB(6),'repair'}));
            % Objective functions
            tempGrid            = GUI.APP(13,[1 6],uigridlayout(obj.app.grid(2),'RowHeight',{'1x'},'ColumnWidth',{135,30,30,'1x'},'Padding',[0 0 0 6],'RowSpacing',0,'ColumnSpacing',5,'BackgroundColor','w'));
            obj.app.titleB(6)   = GUI.APP(1,[1 4],uilabel(tempGrid,'Text',' Objective functions','FontSize',13,'FontColor',[.9 .5 .2],'BackgroundColor',[.9 .9 .9],'FontWeight','bold'));
            obj.app.buttonB(13) = GUI.APP(1,2,uibutton(tempGrid,'Text','Add','FontSize',10,'Tooltip','Add an objective function','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_updateProblem,4}));
            obj.app.buttonB(14) = GUI.APP(1,3,uibutton(tempGrid,'Text','Del','Enable',false,'FontSize',10,'Tooltip','Delete the last objective function','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_updateProblem,-4}));
            obj.pro(1).label    = GUI.APP(14,1,uilabel(obj.app.grid(2),'Text','f1(x) =','HorizontalAlignment','right'));
            obj.pro(1).edit     = GUI.APP(14,[2 4],uieditfield(obj.app.grid(2),'Value','mean(x)','Tooltip','Objective function to be minimized'));
            obj.pro(1).button   = GUI.APP(14,5,uibutton(obj.app.grid(2),'Text','...','Tooltip','Load existing function or data','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_loadFunction,obj.pro(1).edit,{'*.m','MATLAB function';'*.txt;*.dat;*.csv','Text file'}}));
            obj.pro(1).button2  = GUI.APP(14,6,uibutton(obj.app.grid(2),'Text','+','Fontsize',15,'Tooltip','Create new function','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_createFunction,obj.pro(1).edit,'objective'}));
            % Constraint functions
            tempGrid            = GUI.APP(34,[1 6],uigridlayout(obj.app.grid(2),'RowHeight',{'1x'},'ColumnWidth',{135,30,30,'1x'},'Padding',[0 0 0 6],'RowSpacing',0,'ColumnSpacing',5,'BackgroundColor','w'));
            obj.app.titleB(7)   = GUI.APP(1,[1 4],uilabel(tempGrid,'Text',' Constraint functions','FontSize',13,'FontColor',[.9 .5 .2],'BackgroundColor',[.9 .9 .9],'FontWeight','bold'));
            obj.app.buttonB(15) = GUI.APP(1,2,uibutton(tempGrid,'Text','Add','FontSize',10,'Tooltip','Add a constraint function','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_updateProblem,5}));
            obj.app.buttonB(16) = GUI.APP(1,3,uibutton(tempGrid,'Text','Del','FontSize',10,'Tooltip','Delete the last constraint function','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_updateProblem,-5}));
            obj.con(1).labelA   = GUI.APP(35,1,uilabel(obj.app.grid(2),'Text','g1(x) =','HorizontalAlignment','right'));
            obj.con(1).edit     = GUI.APP(35,[2 3],uieditfield(obj.app.grid(2),'Value','0.5-mean(x)','Tooltip','Constraint function to be satisfied'));
            obj.con(1).labelB   = GUI.APP(35,4,uilabel(obj.app.grid(2),'Text','<= 0'));
            obj.con(1).button   = GUI.APP(35,5,uibutton(obj.app.grid(2),'Text','...','Tooltip','Load existing function or data','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_loadFunction,obj.con(1).edit,{'*.m','MATLAB function';'*.txt;*.dat;*.csv','Text file'}}));
            obj.con(1).button2  = GUI.APP(35,6,uibutton(obj.app.grid(2),'Text','+','Fontsize',15,'Tooltip','Create new function','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_createFunction,obj.con(1).edit,'constraint'}));
            % Special difficulties
            tempGrid            = GUI.APP(55,[1 6],uigridlayout(obj.app.grid(2),'RowHeight',{'1x'},'ColumnWidth',{127,30,30,'1x'},'Padding',[0 0 0 6],'RowSpacing',0,'ColumnSpacing',5,'BackgroundColor','w'));
            obj.app.titleB(7)   = GUI.APP(1,[1 4],uilabel(tempGrid,'Text',' Special difficulties','FontSize',13,'FontColor',[.9 .5 .2],'BackgroundColor',[.9 .9 .9],'FontWeight','bold'));
            obj.app.checkB(1)   = GUI.APP(56,[1 6],uicheckbox(obj.app.grid(2),'FontSize',11,'Text','<large> The search space is large, e.g., having more than 100 real decision variables.','WordWrap','on','ValueChangedFcn',@obj.cb_updateFilter));
            obj.app.checkB(2)   = GUI.APP(57,[1 6],uicheckbox(obj.app.grid(2),'FontSize',11,'Text','<expensive> The functions are computationally expensive, e.g., using less than 1000 function evaluations.','WordWrap','on','ValueChangedFcn',@obj.cb_updateFilter));
            obj.app.checkB(3)   = GUI.APP(58,[1 6],uicheckbox(obj.app.grid(2),'FontSize',11,'Text','<multimodal> There are multiple optimal solutions with similar objective values but different decision vectors, all of which should be obtained.','WordWrap','on','ValueChangedFcn',@obj.cb_updateFilter));
            obj.app.checkB(4)   = GUI.APP(59,[1 6],uicheckbox(obj.app.grid(2),'FontSize',11,'Text','<sparse> The optimal solutions are sparse, i.e., most decision variables of which are zero.','WordWrap','on','ValueChangedFcn',@obj.cb_updateFilter));
            
            % The third panel
            obj.app.grid(3)    = GUI.APP([2 3],4,uigridlayout(obj.app.maingrid,'RowHeight',{16,19,16,19,19,16,19,19,19,22,'1.2x','1x'},'ColumnWidth',{'1x','1x','1x'},'Padding',[8 10 8 0],'RowSpacing',3,'ColumnSpacing',5,'BackgroundColor','w'));
            obj.app.labelC(1)  = GUI.APP(1,[1 3],uilabel(obj.app.grid(3),'Text','Number of objectives','VerticalAlignment','bottom','FontSize',12,'FontColor',[.15 .6 .2],'FontWeight','bold'));
            obj.app.stateC(1)  = GUI.APP(2,1,uibutton(obj.app.grid(3),'state','Text','single','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','The problem has a single objective','ValueChangedFcn',{@obj.cb_filter,1}));
            obj.app.stateC(2)  = GUI.APP(2,2,uibutton(obj.app.grid(3),'state','Text','multi','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','The problem has 2 or 3 objectives','ValueChangedFcn',{@obj.cb_filter,2}));
            obj.app.stateC(3)  = GUI.APP(2,3,uibutton(obj.app.grid(3),'state','Text','many','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','The problem has more than 3 objectives','ValueChangedFcn',{@obj.cb_filter,3}));
            obj.app.labelC(2)  = GUI.APP(3,[1 3],uilabel(obj.app.grid(3),'Text','Encoding scheme','VerticalAlignment','bottom','FontSize',12,'FontColor',[.15 .6 .2],'FontWeight','bold'));
            obj.app.stateC(4)  = GUI.APP(4,1,uibutton(obj.app.grid(3),'state','Text','real','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','The decision variables are real numbers','ValueChangedFcn',{@obj.cb_filter,4}));
            obj.app.stateC(5)  = GUI.APP(4,2,uibutton(obj.app.grid(3),'state','Text','integer','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','The decision variables are integers','ValueChangedFcn',{@obj.cb_filter,5}));
            obj.app.stateC(6)  = GUI.APP(4,3,uibutton(obj.app.grid(3),'state','Text','label','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','The decision variables are labels','ValueChangedFcn',{@obj.cb_filter,6}));
            obj.app.stateC(7)  = GUI.APP(5,1,uibutton(obj.app.grid(3),'state','Text','binary','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','The decision variables are binary numbers','ValueChangedFcn',{@obj.cb_filter,7}));
            obj.app.stateC(8)  = GUI.APP(5,2,uibutton(obj.app.grid(3),'state','Text','permutation','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','The decision vector is a permutation','ValueChangedFcn',{@obj.cb_filter,8}));
            obj.app.labelC(3)  = GUI.APP(6,[1 3],uilabel(obj.app.grid(3),'Text','Special difficulties','VerticalAlignment','bottom','FontSize',12,'FontColor',[.15 .6 .2],'FontWeight','bold'));
            obj.app.stateC(9)  = GUI.APP(7,1,uibutton(obj.app.grid(3),'state','Text','large','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','The problem has more than 100 decision variables','ValueChangedFcn',{@obj.cb_filter,9}));
            obj.app.stateC(10) = GUI.APP(7,2,uibutton(obj.app.grid(3),'state','Text','constrained','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','The problem has constraints','ValueChangedFcn',{@obj.cb_filter,10}));
            obj.app.stateC(11) = GUI.APP(7,3,uibutton(obj.app.grid(3),'state','Text','expensive','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','The objectives are computationally time-consuming','ValueChangedFcn',{@obj.cb_filter,11}));
            obj.app.stateC(12) = GUI.APP(8,1,uibutton(obj.app.grid(3),'state','Text','multimodal','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','The objectives are multimodal','ValueChangedFcn',{@obj.cb_filter,12}));
            obj.app.stateC(13) = GUI.APP(8,2,uibutton(obj.app.grid(3),'state','Text','sparse','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','Most decision variables of the optimal solutions are zero','ValueChangedFcn',{@obj.cb_filter,13}));
            obj.app.stateC(14) = GUI.APP(8,3,uibutton(obj.app.grid(3),'state','Text','dynamic','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','The objectives vary periodically','ValueChangedFcn',{@obj.cb_filter,14}));
            obj.app.stateC(15) = GUI.APP(9,1,uibutton(obj.app.grid(3),'state','Text','multitask','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','The problem has multiple tasks to be solved simultaneously','ValueChangedFcn',{@obj.cb_filter,15}));
            obj.app.stateC(16) = GUI.APP(9,2,uibutton(obj.app.grid(3),'state','Text','bilevel','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','The problem has two nested objectives','ValueChangedFcn',{@obj.cb_filter,16}));
            obj.app.stateC(17) = GUI.APP(9,3,uibutton(obj.app.grid(3),'state','Text','robust','FontSize',11,'FontColor',[.15 .6 .2],'BackgroundColor','w','Tooltip','The objectives are influenced by uncertain factors','ValueChangedFcn',{@obj.cb_filter,17}));
            obj.app.labelC(4)  = GUI.APP(10,[1 2],uilabel(obj.app.grid(3),'Text','Algorithms','VerticalAlignment','bottom','FontSize',13,'FontColor',[.2 .4 .7],'FontWeight','bold'));
            obj.app.labelC(5)  = GUI.APP(10,3,uilabel(obj.app.grid(3),'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',10,'FontColor',[.2 .4 .7]));
            obj.app.listC      = GUI.APP(11,[1 3],uilistbox(obj.app.grid(3),'FontColor',[.2 .4 .7],'ValueChangedFcn',{@obj.cb_updateList,3}));
            obj.app.listD      = uilist(obj.app.grid(3),obj.GUI.app.figure,obj.GUI.icon);
            obj.app.tipD       = GUI.APP(12,[1 3],uilabel(obj.app.grid(3),'Text','Select an algorithm','HorizontalAlignment','center'));
            obj.app.listD.grid.Padding = [0,0,0,0];
            GUI.APP(12,[1 3],obj.app.listD.grid);
            
            % The fourth panel
            obj.app.dropD(1) = GUI.APP(1,7,uidropdown(obj.app.maingrid,'BackgroundColor','w','Items',{},'ValueChangedFcn',@obj.cb_slider,'Visible',false));
            obj.app.dropD(2) = GUI.APP(1,7,uidropdown(obj.app.maingrid,'BackgroundColor','w','Items',{},'ValueChangedFcn',@obj.cb_slider));
            obj.app.grid(4)  = GUI.APP([2 3],[6 8],uigridlayout(obj.app.maingrid,'RowHeight',{'1x',40,30},'ColumnWidth',{20,150,'1x','1x',120,30,20},'Padding',[15 10 15 0],'RowSpacing',5,'BackgroundColor','w'));
            obj.app.axes     = GUI.APP(1,[2 6],uiaxes(obj.app.grid(4),'BackgroundColor','w','Box','on'));
            obj.app.waittip  = GUI.APP(1,[2 6],uilabel(obj.app.grid(4),'HorizontalAlignment','center','Text','                 Please wait ... ...','Visible',false));
            tempTb = axtoolbar(obj.app.axes(1),{'rotate','pan','zoomin','zoomout'});
            obj.app.toolD(1)   = axtoolbarbtn(tempTb,'push','Icon',obj.GUI.icon.gif,'Tooltip','Save the evolutionary process to gif','ButtonPushedFcn',@obj.cb_toolbutton1);
            obj.app.toolD(2)   = axtoolbarbtn(tempTb,'push','Icon',obj.GUI.icon.newfigure,'Tooltip','Open in new figure and save to workspace','ButtonPushedFcn',@obj.cb_toolbutton2);
            obj.app.slider     = GUI.APP(2,[1 7],uislider(obj.app.grid(4),'Limits',[0 1],'MajorTicks',0:0.25:1,'MajorTickLabels',{'0%','25%','50%','75%','100%'},'MinorTicks',0:0.01:1,'ValueChangedFcn',@obj.cb_slider));
            obj.app.labelD     = GUI.APP(3,[1 2],uilabel(obj.app.grid(4),'Text','','HorizontalAlignment','left'));
            obj.app.buttonD(1) = GUI.APP(3,3,uibutton(obj.app.grid(4),'push','Text','Start','FontSize',16,'Enable',false,'ButtonpushedFcn',@obj.cb_start));
            obj.app.buttonD(2) = GUI.APP(3,4,uibutton(obj.app.grid(4),'push','Text','Stop','FontSize',16,'Enable',false,'ButtonpushedFcn',@obj.cb_stop));
            obj.app.menuD      = uicontext(obj.GUI.app.figure,120);
            obj.app.menuD.add('  Save best solutions','',{@obj.cb_save,1});
            obj.app.menuD.add('  Save all solutions','',{@obj.cb_save,2});
            obj.app.menuD.flush();
            obj.app.buttonD(3) = GUI.APP(3,[6 7],uibutton(obj.app.grid(4),'push','Text','Save','FontSize',16,'Enable',false,'ButtonpushedFcn',@(~,~)obj.app.menuD.show()));
            
            % Initialization
            obj.cb_updateFilter();
            obj.cb_filter();
        end
    end
    methods(Access = private)
        %% Update the algorithms in the lists
        function cb_filter(obj,~,~,index)
            if nargin > 3
                if index < 4
                    [obj.app.stateC(1:3).Value] = deal(0);
                    obj.app.stateC(index).Value = 1;
                end
            end
            filter = [obj.app.stateC.Value];
            func   = @(s)all(any(repmat([true,filter],size(s,1),1)&s,2)) && all((any(s(:,2:end),1)&filter)==filter);
            % Update the list of algorithms
            show   = cellfun(func,obj.GUI.algList(:,1));
            obj.app.listC.Items = ['(Open File)';obj.GUI.algList(show,2)];
            obj.app.listC.Value = {};
            obj.app.labelC(5).Text = sprintf('%d / %d',sum(show),length(show));
            obj.app.listD.del(1);
            obj.app.tipD.Visible      = 'on';
            obj.app.buttonD(1).Enable = 'off';
            % Update the list of metrics
            show   = cellfun(@(s)func(s(2:end,1:end-2)),obj.GUI.metList(:,1));
            if obj.app.stateC(1).Value == 0 % Multi-objective optimization
                obj.app.dropD(1).Items = ['Population (objectives)';'Population (variables)';obj.GUI.metList(show,2)];
            else                            % Single-objective optimization
                obj.app.dropD(2).Items = ['Population (variables)';obj.GUI.metList(show,2)];
            end
        end
        %% Update the parameter list
        function cb_updateList(obj,~,~,type)
            filename = obj.app.listC.Value;
            if contains(filename,'Open File')
                [file,path] = uigetfile({'*.m','MATLAB class'},'');
                if file ~= 0
                    try
                        filename = fullfile(path,file);
                        f   = fopen(filename);
                        str = fgetl(f);
                        fclose(f);
                        assert(contains(str,'< ALGORITHM'));
                        addpath(path);
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
                case {1,-1}     % Data
                    obj.app.grid(2).RowHeight{8} = max(0,25*sign(index));
                    obj.app.buttonB(1).Enable    = index < 0;
                    obj.app.buttonB(2).Enable    = index > 0;
                    if index > 0
                        str = ',data';
                    else
                        str = '';
                    end
                    obj.app.labelB(7).Text = sprintf('I(N%s) =',str);
                    obj.app.labelB(8).Text = sprintf('R(x%s) =',str);
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
                    obj.app.buttonB(8).Enable     = index < 0;
                    obj.app.buttonB(9).Enable     = index > 0;
                case {4,-4}     % Objective functions
                    if index>0 && length(obj.pro)<20
                        if obj.app.grid(2).RowHeight{8} > 0
                            str = {',data','*data'};
                        else
                            str = {'',''};
                        end
                        item.label   = GUI.APP(length(obj.pro)+14,1,uilabel(obj.app.grid(2),'Text',sprintf('f%d(x%s) =',length(obj.pro)+1,str{1}),'HorizontalAlignment','right'));
                        item.edit    = GUI.APP(length(obj.pro)+14,[2 4],uieditfield(obj.app.grid(2),'Value',sprintf('mean(x%s)',str{2}),'Tooltip','Objective function to be minimized'));
                        item.button  = GUI.APP(length(obj.pro)+14,5,uibutton(obj.app.grid(2),'Text','...','Tooltip','Load existing function or data','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_loadFunction,item.edit,{'*.m','MATLAB function';'*.txt;*.dat;*.csv','Text file'}}));
                        item.button2 = GUI.APP(length(obj.pro)+14,6,uibutton(obj.app.grid(2),'Text','+','Fontsize',15,'Tooltip','Create new function','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_createFunction,item.edit,'objective'}));
                        obj.pro      = [obj.pro,item];
                    elseif index<0 && length(obj.pro)>0
                        delete(obj.pro(end).label);
                        delete(obj.pro(end).edit);
                        delete(obj.pro(end).button);
                        delete(obj.pro(end).button2);
                        obj.pro(end) = [];
                    end
                    obj.app.grid(2).RowHeight(14:33) = num2cell([zeros(1,length(obj.pro))+25,zeros(1,20-length(obj.pro))]);
                    obj.app.buttonB(13).Enable       = length(obj.pro) < 20;
                    obj.app.buttonB(14).Enable       = length(obj.pro) > 1;
                case {5,-5}     % Constraint functions
                    if index>0 && length(obj.con)<20
                        if obj.app.grid(2).RowHeight{8} > 0
                            str = {',data','*data'};
                        else
                            str = {'',''};
                        end
                        item.labelA  = GUI.APP(length(obj.con)+35,1,uilabel(obj.app.grid(2),'Text',sprintf('g%d(x%s) =',length(obj.con)+1,str{1}),'HorizontalAlignment','right'));
                        item.edit    = GUI.APP(length(obj.con)+35,[2 3],uieditfield(obj.app.grid(2),'Value',sprintf('0.5-mean(x%s)',str{2}),'Tooltip','Constraint function to be satisfied'));
                        item.labelB  = GUI.APP(length(obj.con)+35,4,uilabel(obj.app.grid(2),'Text','<= 0'));
                        item.button  = GUI.APP(length(obj.con)+35,5,uibutton(obj.app.grid(2),'Text','...','Tooltip','Load existing function or data','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_loadFunction,item.edit,{'*.m','MATLAB function';'*.txt;*.dat;*.csv','Text file'}}));
                        item.button2 = GUI.APP(length(obj.con)+35,6,uibutton(obj.app.grid(2),'Text','+','Fontsize',15,'Tooltip','Create new function','BackgroundColor','w','ButtonpushedFcn',{@obj.cb_createFunction,item.edit,'constraint'}));
                        obj.con      = [obj.con,item];
                    elseif index<0 && length(obj.con)>0
                        delete(obj.con(end).labelA);
                        delete(obj.con(end).edit);
                        delete(obj.con(end).labelB);
                        delete(obj.con(end).button);
                        delete(obj.con(end).button2);
                        obj.con(end) = [];
                    end
                    obj.app.grid(2).RowHeight(35:54) = num2cell([zeros(1,length(obj.con))+25,zeros(1,20-length(obj.con))]);
                    obj.app.buttonB(15).Enable       = length(obj.con) < 20;
                    obj.app.buttonB(16).Enable       = length(obj.con) > 0;
                case 6          % Evaluation function
                    if strcmp(obj.app.buttonB(10).Text,'Integrate')
                        obj.app.titleB(5).Text      = ' Evaluation function';
                        obj.app.buttonB(8).Enable   = false;
                        obj.app.buttonB(9).Enable   = false;
                        obj.app.buttonB(10).Text    = 'Separate';
                        obj.app.buttonB(10).Tooltip = 'Separate the evaluation function into repair function, objective functions, and constraint functions';
                        obj.app.labelB(8).Text      = '[x,f,g] =';
                        [obj.app.editB(6).Value,obj.app.editB(6).UserData] = deal(obj.app.editB(6).UserData,obj.app.editB(6).Value);
                        obj.app.editB(6).Tooltip = 'Function for evaluating solution';
                        obj.app.grid(2).RowHeight{12}    = 25;
                        obj.app.grid(2).RowHeight(13:54) = {0};
                    else
                        obj.app.titleB(5).Text      = ' Repair function';
                        obj.app.buttonB(9).Enable   = true;
                        obj.app.buttonB(10).Text    = 'Integrate';
                        obj.app.buttonB(10).Tooltip = 'Integrate the repair function, objective functions, and constraint functions into a single evaluation function';
                        if obj.app.grid(2).RowHeight{8} > 0
                            str = ',data';
                        else
                            str = '';
                        end
                        obj.app.labelB(8).Text = sprintf('R(x%s) =',str);
                        [obj.app.editB(6).Value,obj.app.editB(6).UserData] = deal(obj.app.editB(6).UserData,obj.app.editB(6).Value);
                        obj.app.editB(6).Tooltip = 'Function for repairing invalid solution';
                        obj.app.grid(2).RowHeight(13:13+length(obj.pro)) = {25};
                        obj.app.grid(2).RowHeight(34:34+length(obj.con)) = {25};
                    end
            end
            if nargin < 5 || ~noFlush
                obj.cb_updateFilter();
            end
        end
        %% Validate the problem
        function PRO = cb_validation(obj,para,paraname)
            if nargin < 2
                para     = {1,1};
                paraname = 'maxFE';
            end
            Str = GetSetProblem(obj);
            try
                if length(Str) > 6
                    PRO = UserProblem('N',para{1},paraname,para{2},'encoding',Str{1},'lower',Str{2},'upper',Str{3},'data',Str{4},'initFcn',Str{5},'decFcn',Str{6},'objFcn',Str{7},'conFcn',Str{8});
                else
                    PRO = UserProblem('N',para{1},paraname,para{2},'encoding',Str{1},'lower',Str{2},'upper',Str{3},'data',Str{4},'initFcn',Str{5},'evalFcn',Str{6});
                end
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
                        Str = {'ones(1,10)','zeros(1,10)','ones(1,10)','','','','mean(x.^2)',''};
                    case 3
                        Str = {'ones(1,15)','zeros(1,15)','ones(1,15)','','','',{'x(1)+mean(x(2:end).^2)','1-x(1).^2+mean(x(2:end).^2)'},''};
                    case 4
                        Str = {'ones(1,20)','zeros(1,20)','ones(1,20)','','','',{'x(1)*x(2)*x(3)+mean(x(4:end).^2)','x(1)*x(2)*(1-x(3))+mean(x(4:end).^2)','x(1)*(1-x(2))+mean(x(4:end).^2)','1-x(1)+mean(x(4:end).^2)'},''};
                    case 5
                        Str = {'ones(1,25)','zeros(1,25)','ones(1,25)','','','','mean(x.^2)','0.2-mean(x.^2)'};
                    case 6
                        Str = {'ones(1,30)','zeros(1,30)','ones(1,30)','','','',{'x(1)+mean(x(2:end).^2)','1-x(1).^2+mean(x(2:end).^2)'},{'0.2-mean(x(2:end).^2)','mean(x(2:end).^2)-0.8'}};
                    case 7
                        Str = {'ones(1,35)','zeros(1,35)','ones(1,35)','rand(35)','','','mean((x*data).^2)',''};
                    case 8
                        Str = {'zeros(1,40)+2','zeros(1,40)','zeros(1,40)+10','','','','mean((x-5).^2)',''};
                    case 9
                        Str = {'zeros(1,100)+4','','','','rand(N,100)<0.2','','1-mean(x)','mean(x)-0.1'};
                    case 10
                        Str = {'zeros(1,20)+5','','','','','','mean(cumsum(movsum(x,3)))',''};
                    case 11
                        Str = {'1,1,1,2,2,2,4,4,4','zeros(1,9)','zeros(1,9)+10','','','','sum((x-pi).^2)',''};
                end
                if isempty(obj.app.dropA.UserData)
                    obj.app.dropA.UserData = GetSetProblem(obj);
                end
                GetSetProblem(obj,Str);
            end
        end
        %% Load a problem
        function cb_loadProblem(obj,~,~)
            [file,folder] = uigetfile({'*.m','MATLAB code'});
            if ischar(file)
                try
                    f = fopen(fullfile(folder,file));
                    fgetl(f); fgetl(f); fgetl(f); fgetl(f); fgetl(f); fgetl(f);
                    Str{1} = regexp(fgetl(f),'(?<=''encoding'','').*(?='',...)','match','once');
                    Str{2} = regexp(fgetl(f),'(?<=''lower'','').*(?='',...)','match','once');
                    Str{3} = regexp(fgetl(f),'(?<=''upper'','').*(?='',...)','match','once');
                    Str{4} = regexp(fgetl(f),'(?<=''data'','').*(?='',...)','match','once');
                    Str{5} = regexp(fgetl(f),'(?<=''initFcn'','').*(?='',...)','match','once');
                    str    = fgetl(f);
                    if contains(str,'decFcn')
                        Str{6} = regexp(str,'(?<=''decFcn'','').*(?='',...)','match','once');
                        Str{7} = eval(regexp(fgetl(f),'(?<=''objFcn'',).*(?=,...)','match','once'));
                        Str{8} = eval(regexp(fgetl(f),'(?<=''conFcn'',).*(?=\);)','match','once'));
                    else
                        Str{6} = regexp(str,'(?<=''evalFcn'','').*(?=''\);)','match','once');
                    end
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
            [Name,Path] = uiputfile({'*.m','MATLAB code'},'','new');
            if ischar(Name)
                try
                    Str      = GetSetProblem(obj);
                    [~,name] = fileparts(Name);
                    Code = {['classdef ',name,' < UserProblem % < PROBLEM']
                           '% The code is automatically generated and should not be modified'
                           ''
                           '    methods'
                           ['        function obj = ',name,'(varargin)']
                           '            obj = obj@UserProblem(varargin{:},...'
                           ['            ''encoding'',''',Str{1},''',...']
                           ['            ''lower'',''',Str{2},''',...']
                           ['            ''upper'',''',Str{3},''',...']
                           ['            ''data'',''',Str{4},''',...']
                           ['            ''initFcn'',''',Str{5},''',...']};
                    if length(Str) > 6
                        Code = [Code;
                               ['            ''decFcn'',''',Str{6},''',...']
                               ['            ''objFcn'',{',strjoin(strcat('''',Str{7},''''),','),'},...']
                               ['            ''conFcn'',{',strjoin(strcat('''',Str{8},''''),','),'});']];
                    else
                        Code = [Code;
                               ['            ''evalFcn'',''',Str{6},''');']];
                    end
                    Code = [Code;
                            '        end'
                            '    end'
                            'end'];
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
        %% Load existing data or function
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
        %% Create new function
        function cb_createFunction(obj,~,~,h,type)
            if length(h.Value)>2 && h.Value(1)=='<' && h.Value(end)=='>'
                [file,folder] = uiputfile({'*.m','MATLAB function'},'',fullfile(fileparts(h.Value(2:end-1)),'myFcn'));
            else
                [file,folder] = uiputfile({'*.m','MATLAB function'},'',fullfile(cd,'myFcn'));
            end
            if ischar(file)
                try
                    fid      = fopen(fullfile(folder,file),'wt');
                    [~,name] = fileparts(file);
                    switch type
                        case 'initialization'           % Initialization function
                            if obj.app.grid(2).RowHeight{8} > 0
                                str = ',data';
                            else
                                str = '';
                            end
                            fprintf(fid,'function PopDec = %s(N%s)\n',name,str);
                            fprintf(fid,'%% This is a user-defined %s function\n',type);
                            fprintf(fid,'%% Input ''N'' is a scalar denoting the number of generated solutions\n');
                            if obj.app.grid(2).RowHeight{8} > 0
                                fprintf(fid,'%% Input ''data'' can be any type of user-defined data\n');
                            end
                            fprintf(fid,'%% Output ''PopDec'' is a matrix with each row denoting the decision vector of a solution\n\n');
                            fprintf(fid,'    PopDec = rand(N,10);\n');
                            fprintf(fid,'end');
                        case 'repair'                   % Repair function or evaluation function
                            if obj.app.grid(2).RowHeight{8} > 0
                                str = {',data','*data'};
                            else
                                str = {'',''};
                            end
                            if strcmp(obj.app.buttonB(10).Text,'Integrate')
                                fprintf(fid,'function x = %s(x%s)\n',name,str{1});
                            else
                                fprintf(fid,'function [x,f,g] = %s(x%s)\n',name,str{1});
                                type = 'evaluation';
                            end
                            fprintf(fid,'%% This is a user-defined %s function\n',type);
                            fprintf(fid,'%% Input ''x'' is a row vector denoting the decision vector of a solution\n');
                            if obj.app.grid(2).RowHeight{8} > 0
                                fprintf(fid,'%% Input ''data'' can be any type of user-defined data\n');
                            end
                            if strcmp(obj.app.buttonB(10).Text,'Integrate')
                                fprintf(fid,'%% Output ''x'' denotes the repaired decision vector of the solution\n\n');
                            else
                                fprintf(fid,'%% Output ''x'' denotes the repaired decision vector of the solution\n');
                                fprintf(fid,'%% Output ''f'' denotes all the objective values of the solution\n');
                                fprintf(fid,'%% Output ''g'' denotes all the constraint values of the solution\n\n');
                            end
                            if strcmp(obj.app.buttonB(10).Text,'Integrate')
                                fprintf(fid,'    x = max(min(x,1),0);\n');
                            else
                                fprintf(fid,'    x    = max(min(x,1),0);\n');
                                fprintf(fid,'    f(1) = mean(x%s);\n',str{2});
                                fprintf(fid,'    f(2) = 1-mean(x);\n');
                                fprintf(fid,'    g(1) = mean(x%s)-0.8;\n',str{2});
                                fprintf(fid,'    g(2) = 0.2-mean(x%s);\n',str{2});
                            end
                            fprintf(fid,'end');
                        case {'objective','constraint'} % Objective and constraint functions
                            if obj.app.grid(2).RowHeight{8} > 0
                                str = {',data','*data'};
                            else
                                str = {'',''};
                            end
                            fprintf(fid,'function y = %s(x%s)\n',name,str{1});
                            fprintf(fid,'%% This is a user-defined %s function\n',type);
                            fprintf(fid,'%% Input ''x'' is a row vector denoting the decision vector of a solution\n');
                            if obj.app.grid(2).RowHeight{8} > 0
                                fprintf(fid,'%% Input ''data'' can be any type of user-defined data\n');
                            end
                            fprintf(fid,'%% Output ''y'' is a scalar denoting the %s value of the solution\n',type);
                            if strcmp(type,'objective')
                                fprintf(fid,'%% The objective is to be minimized\n\n');
                                fprintf(fid,'    y = mean(x%s);\n',str{2});
                            else
                                fprintf(fid,'%% The constraint is satisfied if and only it is not positive\n\n');
                                fprintf(fid,'    y = mean(x%s)-0.5;           %% Constraint mean(x%s)<=0.5\n',str{2},str{2});
                                fprintf(fid,'    y = 0.5-mean(x%s);           %% Constraint mean(x%s)>=0.5\n',str{2},str{2});
                                fprintf(fid,'    y = abs(0.5-mean(x%s))-1e-6; %% Constraint mean(x%s)==0.5\n',str{2},str{2});
                            end
                            fprintf(fid,'end');
                    end
                    fclose(fid);
                    web(['file://',fullfile(folder,file)],'-browser');
                    h.Value = sprintf('<%s>',fullfile(folder,file));
                catch err
                    uialert(obj.GUI.app.figure,'Fail to create a new function, please refer to the command window for details.','Error');
                    rethrow(err);
                end
            end
        end
        %% Update the filter
        function cb_updateFilter(obj,~,~)
            oldState = [obj.app.stateC.Value];
            if strcmp(obj.app.buttonB(10).Text,'Integrate')
                obj.app.stateC(1).Value = length(obj.pro) <= 1;
                obj.app.stateC(2).Value = length(obj.pro) >= 2 && length(obj.pro) <= 3;
                obj.app.stateC(3).Value = length(obj.pro) >= 4;
            end
            encoding = str2num(obj.app.editB(1).Value);
            obj.app.stateC(4).Value = any(encoding==1);
            obj.app.stateC(5).Value = any(encoding==2);
            obj.app.stateC(6).Value = any(encoding==3);
            obj.app.stateC(7).Value = any(encoding==4);
            obj.app.stateC(8).Value = any(encoding==5);
            obj.app.stateC(9).Value = obj.app.checkB(1).Value;
            if strcmp(obj.app.buttonB(10).Text,'Integrate')
                obj.app.stateC(10).Value = ~isempty(obj.con);
            end
            obj.app.stateC(11).Value = obj.app.checkB(2).Value;
            obj.app.stateC(12).Value = obj.app.checkB(3).Value;
            obj.app.stateC(13).Value = obj.app.checkB(4).Value;
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
                PRO = obj.cb_validation(para(1:2),item.label(2).Text);
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
                [obj.app.stateC.Enable]       = deal('off');
                obj.app.listC.Enable          = 'off';
                obj.app.listD.Enable          = 'off';
                [obj.app.toolD(1:2).Visible]  = deal('off');
                obj.app.buttonD(1).Text       = 'Pause';
                obj.app.buttonD(2).Enable     = 'on';
                obj.app.buttonD(3).Enable     = 'on';
                if PRO.M > 1
                    obj.app.dropD(1).Value   = obj.app.dropD(1).Items{1};
                    obj.app.dropD(1).Visible = 'on';
                    obj.app.dropD(2).Visible = 'off';
                else
                    obj.app.dropD(2).Value   = obj.app.dropD(2).Items{1};
                    obj.app.dropD(1).Visible = 'off';
                    obj.app.dropD(2).Visible = 'on';
                end
                obj.app.dropD(1).Enable = 'off';
                obj.app.dropD(2).Enable = 'off';
                set([obj.pro.edit,obj.con.edit],'Enable','off');
                set([obj.pro.button,obj.pro.button2,obj.con.button,obj.con.button2],'Enable','off');
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
            [obj.app.stateC.Enable]       = deal('on');
            obj.app.listC.Enable          = 'on';
            obj.app.listD.Enable          = 'on';
            [obj.app.toolD(1:2).Visible]  = deal('on');
            obj.app.buttonD(1).Text       = 'Start';
            obj.app.buttonD(2).Enable     = 'off';
            obj.app.dropD(1).Enable       = 'on';
            obj.app.dropD(2).Enable       = 'on';
            set([obj.pro.edit,obj.con.edit],'Enable','on');
            set([obj.pro.button,obj.pro.button2,obj.con.button,obj.con.button2],'Enable','on');
        end
        %% Save the result
        function cb_save(obj,~,~,type)
            ALG   = obj.data{1};
            PRO   = obj.data{2};
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
                    switch obj.app.dropD(1).Value
                        case 'Population (objectives)'
                            PRO.DrawObj(ALG.result{index,2});
                        case 'Population (variables)'
                            PRO.DrawDec(ALG.result{index,2});
                        otherwise
                            obj.app.waittip.Visible = 'on'; drawnow();
                            Draw(ALG.CalMetric(obj.app.dropD(1).Value),'-k.','LineWidth',1.5,'MarkerSize',10,{'Number of function evaluations',strrep(obj.app.dropD(1).Value,'_',' '),[]});
                            obj.app.waittip.Visible = 'off';
                    end
                else
                    % Show the result of single-objective optimization
                    switch obj.app.dropD(2).Value
                        case 'Population (variables)'
                            PRO.DrawDec(ALG.result{index,2});
                        otherwise
                            Draw(ALG.CalMetric(obj.app.dropD(2).Value),'-k.','LineWidth',1.5,'MarkerSize',10,{'Number of function evaluations',strrep(obj.app.dropD(2).Value,'_',' '),[]});
                    end
                end
            end
        end
        %% Create the gif
        function cb_toolbutton1(obj,~,~)
            if ~isempty(obj.data)
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
    end
end

function Str = GetSetProblem(obj,Str)
    if nargout > 0  % Get the problem
        Str{1} = obj.app.editB(1).Value;              	% Encoding
        Str{2} = obj.app.editB(2).Value;             	% Lower bound
        Str{3} = obj.app.editB(3).Value;              	% Upper bound
        if obj.app.grid(2).RowHeight{8} > 0          	% Data
            Str{4} = obj.app.editB(4).Value;
        else
            Str{4} = '';
        end
        if obj.app.grid(2).RowHeight{10} > 0          	% Initialization function
            Str{5} = obj.app.editB(5).Value;
        else
            Str{5} = '';
        end
        if obj.app.grid(2).RowHeight{12} > 0           	% Repair function or evaluation function
            Str{6} = obj.app.editB(6).Value;
        else
            Str{6} = '';
        end
        if strcmp(obj.app.buttonB(10).Text,'Integrate')
            Str{7} = get([obj.pro.edit],'Value');     	% Objective functions
            if ~iscell(Str{7}); Str{7} = Str(7); end
            Str{8} = get([obj.con.edit],'Value');     	% Constraint functions
            if ~iscell(Str{8}); Str{8} = Str(8); end
        end
    end
    if nargin > 1	% Set the problem
        obj.app.editB(1).Value  = Str{1};               % Number of decision variables
        obj.app.editB(2).Value  = Str{2};               % Lower bound
        obj.app.editB(3).Value  = Str{3};            	% Upper bound
        obj.app.editB(4).Value  = Str{4};               % Data
        obj.cb_updateProblem([],[],1*(-1)^isempty(obj.app.editB(4).Value),1);
        obj.app.editB(5).Value  = Str{5};            	% Initialization function
        obj.cb_updateProblem([],[],2*(-1)^isempty(obj.app.editB(5).Value),1);
        if length(Str) > 6
            if ~strcmp(obj.app.buttonB(10).Text,'Integrate')
                obj.cb_updateProblem([],[],6,1);
            end
            obj.app.editB(6).Value = Str{6};         	% Repair function
            obj.cb_updateProblem([],[],3*(-1)^isempty(obj.app.editB(6).Value),1);
            for i = 1 : length(obj.pro)                	% Objective functions
                obj.cb_updateProblem([],[],-4,1);
            end
            if ~iscell(Str{7}); Str{7} = Str(7); end
            Str{7}(cellfun(@isempty,Str{7})) = [];
            for i = 1 : length(Str{7})
                obj.cb_updateProblem([],[],4,1);
                obj.pro(end).edit.Value = Str{7}{i};
            end
            for i = 1 : length(obj.con)             	% Constraint functions
                obj.cb_updateProblem([],[],-5,1);
            end
            if ~iscell(Str{8}); Str{8} = Str(8); end
            Str{8}(cellfun(@isempty,Str{8})) = [];
            for i = 1 : length(Str{8})
                obj.cb_updateProblem([],[],5,1);
                obj.con(end).edit.Value = Str{8}{i};
            end
        else
            if strcmp(obj.app.buttonB(10).Text,'Integrate')
                obj.cb_updateProblem([],[],6,1);
            end
            obj.app.editB(6).Value = Str{6};            % Evaluation function
        end
        obj.cb_updateFilter();
    end
end