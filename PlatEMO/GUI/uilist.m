classdef uilist < handle
%uilist - Panel of parameter setting.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        Enable;     % State of the objects
    end
    properties(SetAccess = private)
        grid;       % The main panel
        items;      % The setting groups
        menus;      % The menu objects
        current;    % Index of current item
    end
    methods
        %% Constructor
        function obj = uilist(parent,fig,icon)
            obj.grid  = uigridlayout(parent,'RowHeight',{20},'ColumnWidth',{15,'1x','1x',20},'Padding',[5 0 5 0],'RowSpacing',5,'ColumnSpacing',3,'Scrollable','on','BackgroundColor','w');
            obj.menus = uicontext(fig,100);
            obj.menus.add('Open file',icon.file,@obj.cb_openfile);
            obj.menus.add('Open folder',icon.folder,@obj.cb_openfolder);
            obj.menus.add('Search online',icon.scholar,@obj.cb_search);
            obj.menus.flush();
            obj.menus(2) = uicontext(fig,100);
            obj.menus(2).add('Move up',icon.moveup,@obj.cb_moveup);
            obj.menus(2).add('Move down',icon.movedown,@obj.cb_movedown);
            obj.menus(2).add('Delete',icon.delete,@obj.cb_delete);
            obj.menus(2).add('Open file',icon.file,@obj.cb_openfile);
            obj.menus(2).add('Open folder',icon.folder,@obj.cb_openfolder);
            obj.menus(2).add('Search online',icon.scholar,@obj.cb_search);
            obj.menus(2).gaps(3).Visible = false;
            obj.menus(2).flush();
        end
        %% Add a new item
        function add(obj,filename,type)
            % Read the comments in the head of the .m file
            strSet = {};
            if isempty(which(filename))
                return;
            end
            [~,name] = fileparts(filename);
            if ~isempty(obj.items) && ismember(name,get([obj.items.title],{'Text'}))
                return;
            end
            f = fopen(filename);
            fgetl(f);
            singleobj = contains(fgetl(f),'<single>');
            % Read the summary and parameter settings
            str = fgetl(f);
            while ischar(str) && ~isempty(regexp(str,'^\s*%\s*','once'))
                strSet = [strSet,{regexprep(str,'^\s*%\s*','','once')}];
                str    = fgetl(f);
            end
            % Read the reference
            while ischar(str) && isempty(regexp(str,'---\s*Reference\s*---','once'))
                str = fgetl(f);
            end
            str  = fgetl(f);
            cite = [];
            while ischar(str) && isempty(regexp(str,'^\s*%\s*---','once'))
                cite = [cite,' ',regexprep(str,'^\s*%\s*','','once')];
                str  = fgetl(f);
            end
            fclose(f); 
            % Obtain all the parameter settings
            loc = 1;
            while loc <= length(strSet) && isempty(regexp(strSet{loc},'\s*---\s*','once'))
                loc = loc + 1;
            end
            Comment = strjoin(strSet(1:loc-1),' ');
            if type == 3        % For algorithms in application module
                Parameter = {'N','100','Population size';'maxFE','10000','Maximum number of function evaluations'};
                type = 1;
            elseif type > 0     % For algorithms in other modules
                Parameter = {};
            else                % For problems in other modules
                Parameter = {'N','100','Population size'
                             'M','','Number of objectives'
                             'D','','Number of decision variables'
                             'maxFE','10000','Maximum number of function evaluations'};
            end
            for i = loc : length(strSet)
                str = regexp(strSet{i},'\s*---\s*','split');
                if length(str) >= 2
                    if length(str) == 2
                        str{3} = '';
                    end
                    Parameter = [Parameter;str(1:3)];
                end
            end
            % Generate all the items
            if type > 0
                color = [.2 .4 .7];
            else
                color = [.9 .5 .2];
            end
            item.cite  = cite;
            item.type  = type;
            item.fold  = false;
            item.title = uibutton(obj.grid,'Text',name,'FontSize',11,'BackgroundColor',color,'FontColor',[1 1 1],'Tooltip',Comment,'ButtonPushedFcn',@obj.cb_fold);
            item.title.Layout.Column = [1 4];
            item.tip = uibutton(obj.grid,'Text','¨‹','FontSize',10,'BackgroundColor',color,'FontColor',[1 1 1],'Tooltip',Comment,'ButtonPushedFcn',@obj.cb_callmenu);
            item.tip.Layout.Column = 4;
            item.label  = [];
            item.edit   = [];
            item.button = [];
            for i = 1 : size(Parameter,1)
                item.label = [item.label,uilabel(obj.grid,'Text',Parameter{i,1},'Tooltip',Parameter{i,3})];
                item.label(end).Layout.Column = 2;
                if type<0 && i==2 && singleobj  % Single-objective optimization
                    item.edit = [item.edit,uieditfield(obj.grid,'Value','1','Enable','off','HorizontalAlignment','right','Tooltip',Parameter{i,3},'UserData',false)];
                else                            % Multi-objective optimization
                    item.edit = [item.edit,uieditfield(obj.grid,'Value',Parameter{i,2},'HorizontalAlignment','right','Tooltip',Parameter{i,3},'UserData',true)];
                end
                item.edit(end).Layout.Column = 3;
                if strcmp(Parameter{i,1},'maxFE')
                    item.button = uibutton(obj.grid,'Text','¨‹','FontSize',10,'BackgroundColor','w','UserData',{item.label(end),item.edit(end)},'ButtonPushedFcn',@obj.cb_switch);
                    item.button.Layout.Column = 4;
                end
            end
            obj.items = [obj.items,item];
        end
        %% Delete existing items
        function del(obj,index,type)
            if ~isempty(obj.items)
                if nargin >= 3
                    index = find([obj.items.type]==type);
                end
                for i = index
                    delete(obj.items(i).title);
                    delete(obj.items(i).tip);
                    delete([obj.items(i).label]);
                    delete([obj.items(i).edit]);
                    delete(obj.items(i).button);
                end
                obj.items(index) = [];
            end
        end
        %% Flush the list
        function flush(obj)
            if ~isempty(obj.items)
                % Sort the items
                [~,rank]  = sort([obj.items.type],'descend');
                obj.items = obj.items(rank);
                % Relocate the items
                obj.grid.RowHeight = repmat({22},1,length(obj.items)+sum(arrayfun(@(s)length(s.label),obj.items).*~[obj.items.fold]));
                loc = 0;
                for i = 1 : length(obj.items)
                    loc = loc + 1;
                    obj.items(i).title.Layout.Row = loc;
                    obj.items(i).tip.Layout.Row   = loc;
                    if ~obj.items(i).fold
                        for j = 1 : length(obj.items(i).label)
                            loc = loc + 1;
                            obj.items(i).label(j).Layout.Row = loc;
                            obj.items(i).edit(j).Layout.Row  = loc;
                        end
                        if ~isempty(obj.items(i).button)
                            obj.items(i).button.Layout.Row = obj.items(i).button.UserData{1}.Layout.Row;
                        end
                    end
                end
            end
        end
        %% Change the state of the objects
        function set.Enable(obj,value)
            set([obj.items.title],'Enable',value);
            set([obj.items.tip],'Enable',value);
            set([obj.items.button],'Enable',value);
            edits = [obj.items.edit];
            set(edits,'Enable',value);
            set(edits(~cell2mat(get(edits,'UserData'))),'Enable','off');
        end
    end
    methods(Access = private)
        %% Fold or unfold the item
        function cb_fold(obj,ui,~)
            for i = 1 : length(obj.items)
                if obj.items(i).title == ui
                    obj.items(i).fold = ~obj.items(i).fold;
                    if obj.items(i).fold
                        obj.items(i).tip.Text = '¡ø';
                    else
                        obj.items(i).tip.Text = '¨‹';
                    end
                    if ~isempty(obj.items(i).label)
                        [obj.items(i).label.Visible] = deal(~obj.items(i).fold);
                        [obj.items(i).edit.Visible]  = deal(~obj.items(i).fold);
                        if ~isempty(obj.items(i).button)
                            obj.items(i).button.Visible = ~obj.items(i).fold;
                        end
                        obj.flush();
                    end
                    break;
                end
            end
        end
        %% Show the menu
        function cb_callmenu(obj,ui,~)
            for i = 1 : length(obj.items)
                if obj.items(i).tip == ui
                    obj.current = i;
                    obj.menus(abs(obj.items(i).type)).show();
                    break;
                end
            end
        end
        %% Switch between maxFE and maxRuntime
        function cb_switch(obj,ui,~)
            if strcmp(ui.UserData{1}.Text,'maxFE')
                set(ui.UserData{1},'Text','maxRuntime','Tooltip','Maximum runtime (in second)');
                set(ui.UserData{2},'Value','1','Tooltip','Maximum runtime (in second)');
                set(ui,'Text','¡ø');
            else
                set(ui.UserData{1},'Text','maxFE','Tooltip','Maximum number of function evaluations');
                set(ui.UserData{2},'Value','10000','Tooltip','Maximum number of function evaluations');
                set(ui,'Text','¨‹');
            end
        end
        %% Open file
        function cb_openfile(obj,ui,~)
            ui.Parent.Visible = false;
            web(['file://',which(obj.items(obj.current).title.Text)],'-browser');
        end
        %% Open folder
        function cb_openfolder(obj,ui,~)
            ui.Parent.Visible = false;
            web(['file://',fileparts(which(obj.items(obj.current).title.Text))],'-browser');
        end
        %% Search online
        function cb_search(obj,ui,~)
            ui.Parent.Visible = false;
            web(['https://scholar.google.com/scholar?q=%',strjoin(cellstr(dec2hex(double(obj.items(obj.current).cite))),'%')],'-browser');
        end
        %% Move the item up
        function cb_moveup(obj,ui,~)
            ui.Parent.Visible = false;
            if obj.current > 1
                obj.items = obj.items([1:obj.current-2,obj.current,obj.current-1,obj.current+1:end]);
                obj.flush();
            end
        end
        %% Move the item down
        function cb_movedown(obj,ui,~)
            ui.Parent.Visible = false;
            if obj.current < length(obj.items)
                obj.items = obj.items([1:obj.current-1,obj.current+1,obj.current,obj.current+2:end]);
                obj.flush();
            end
        end
        %% Delete the item
        function cb_delete(obj,ui,~)
            ui.Parent.Visible = false;
            obj.del(obj.current);
            obj.flush();
        end
    end
end