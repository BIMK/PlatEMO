classdef ParameterList_Item < newPanel
%ParameterList_Item - The object used in ParameterList.
%
%   This is the class of item used in ParameterList, which cannot be
%   instantiated independently.
%
%   See also ParameterList

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = ?newGUI)
        type;           % Type of menu (0:none, 1:menu, 2:menu with delete option, negative: including parameters of environment)
        title;          % The title button
        tip;        	% The menu button
        labels;         % All the labels
        edits;          % All the edit boxes
        name;           % Name of the function
        summary = [];	% Summary of the function
        ref     = [];   % Reference of the function
        show = true;    % Whether show the labels and edits
    end
    methods(Access = ?ParameterList)
        %% Constructor
        function obj = ParameterList_Item(parent,mfileName,color,type)
            obj@newPanel(parent,[1 1 parent.position(3) 1,1 0 1 0],parent.color);
            set(obj,'name',mfileName,'type',type);
            obj.createObj(obj.strSet2PStruct(obj.readMfile(mfileName)),color);
        end
        %% Read the comments in the head of the .m file
        function strSet = readMfile(obj,mfileName)
            strSet    = {};
            mfileName = [mfileName,'.m'];
            if ~isempty(which(mfileName))
                f = fopen(mfileName);
                fgetl(f);
                fgetl(f);
                % Read the summary and parameter settings
                str = fgetl(f);
                while ischar(str) && ~isempty(regexp(str,'^\s*%\s*','once'))
                    strSet = [strSet,{regexprep(str,'^\s*%\s*','','once')}];
                    str    = fgetl(f);
                end
                % Read the reference
                while ischar(str) && isempty(strfind(str,'--- Reference ---'))
                    str = fgetl(f);
                end
                str = fgetl(f);
                while ischar(str) && isempty(strfind(str,'--- Copyright ---'))
                    obj.ref = [obj.ref,' ',regexprep(str,'^\s*%\s*','','once')];
                    str     = fgetl(f);
                end
                fclose(f);
            end
        end
        %% Convert the comments to parameter struct
        function Parameter = strSet2PStruct(obj,strSet)
            % Generate the parameter struct
            current = 1;
            while current <=length(strSet) && isempty(regexp(strSet{current},'\s*---\s*','once'))
                current = current + 1;
            end
            obj.summary = strjoin(strSet(1:current-1),' ');
            if obj.type > 0
                Parameter = {};
            else
                Parameter = {'N','100','Population size'
                             'M','','Number of objectives'
                             'D','','Number of decision variables'
                             'evaluation','10000','Number of evaluations'};
            end
            for i = current : length(strSet)
                str = regexp(strSet{i},'\s*---\s*','split');
                if length(str) >= 2
                    if length(str) == 2
                        str{3} = '';
                    end
                    Parameter = [Parameter;str(1:3)];
                end
            end
        end
        %% Create all the GUI objects
        function createObj(obj,Parameter,color)
            obj.position(4) = (1+size(Parameter,1))*25;
            rate  = obj.position(3)/180;
            start = obj.position(4) + 8;
            % Create the items
            for i = 1 : size(Parameter,1)
                obj.labels = [obj.labels,newLabel(obj,[5*rate start-25*(i+1) 100*rate 17,0 0 1 0],Parameter{i,1},'HorizontalAlignment','left','FontSize',9,'movecallback',@(~,~)obj.cb_move(Parameter{i,3}))];
                obj.edits  = [obj.edits,newEdit(obj,[100*rate start-25*(i+1) 75*rate 18,0 0 1 0],'String',Parameter{i,2},'FontSize',9,'movecallback',@(~,~)obj.cb_move(Parameter{i,3}))];
            end
            % Create the title button
            obj.title = newLabel(obj,[1 start-25 obj.position(3) 18,1 1 1 0],obj.name,'BackgroundColor',color,'ForegroundColor',[1 1 1],'callback',@obj.cb_button,'movecallback',@(~,~)obj.cb_move(obj.summary));
            % Create the menu button
            if obj.type ~= 0
                obj.tip = newButtonSpecial(obj,[165*rate start-25 12*rate 18,0 0 1 0],1,color,'ForegroundColor',[1 1 1],'callback',@obj.cb_tip,'movecallback',@(~,~)obj.cb_move(obj.summary));
            end
        end
    end
    methods(Access = protected)
        %% The callback of obj
        function cb_button(obj,hObject,eventdata)
            obj.show = ~obj.show;
            if ~obj.show
                obj.position(4)       = 19;
                obj.title.position(2) = 2;
                if ~isempty(obj.tip)
                    obj.tip.position(2) = 2;
                end
                if ~isempty(obj.labels)
                    [obj.labels.visible] = deal(false);
                    [obj.edits.visible]  = deal(false);
                end
            else
                obj.position(4)       = (1+length(obj.labels))*25;
                obj.title.position(2) = obj.position(4)+8-25;
                if ~isempty(obj.tip)
                    obj.tip.position(2) = obj.position(4)+8-25;
                end
                if ~isempty(obj.labels)
                    [obj.labels.visible] = deal(true);
                    [obj.edits.visible]  = deal(true);
                end
            end
            obj.parent.cb_updateList();
        end
        %% The callback of obj.tip
        function cb_tip(obj,hObject,eventdata)
            obj.parent.currentIndex = find([obj.parent.items.handle]==obj.handle,1);
            if abs(obj.type) == 2
                obj.parent.menus(2).items(1).enable = obj.parent.currentIndex > 1;
                obj.parent.menus(2).items(2).enable = obj.parent.currentIndex < length(obj.parent.items);
            end
            obj.parent.menus(abs(obj.type)).show();
        end
        %% The movecallback of all the other controls in obj
        function cb_move(obj,str)
            obj.parent.note.handle.String = str;
        end
    end
end