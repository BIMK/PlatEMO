classdef GLOBAL < handle
%GLOBAL - The class of the experimental setting
%
%   This is the class of the experimental setting. An object of GLOBAL
%   class stores all the propeties used in the algorithm, including the
%   population size, number of objectives, maximum number of evaluations
%   and so on. The object also has several methods which can be invoked by
%   algorithm, problem and operator functions.
%
% GLOBAL properties:
%   N               <public>	population size
%   M               <read-only>	number of objectives
%   D               <read-only>	number of variables
%   lower           <read-only>	lower bound of each decision variable
%   upper           <read-only>	upper bound of each decision variable
%   algorithm       <read-only>	algorithm function
%   problem         <read-only>	problem function
%   operator        <read-only>	operator function
%   evaluated       <read-only>	number of evaluated individuals
%   evaluation      <read-only>	maximum number of evaluations
%   gen             <read-only>	current generation
%   maxgen          <read-only>	maximum generation
%   run             <read-only>	run No.
%   runtime         <read-only>	runtime
%   result          <read-only>	set of recorded populations
%   PF              <read-only>	true Pareto front
%   parameter       <read-only>	specified parameters for algorithm, problem and operator
%   mode            <read-only>	run mode (1.show result 2.save result 3.run outputFcn)
%   outputFcn    	<read-only>	function invoked after each generation when mode = 3
%
% GLOBAL methods:
%   GLOBAL          <public>    the constructor, all the properties will be
%                               set when the object is creating
%   Start           <public>    start running the algorithm
%   Initialization  <public>    generate an initial random population
%   NotTermination  <public>	terminate the running if the number of
%                               evaluations has exceeded
%   Variation       <public>    generate offspring population
%   VariationDec    <public>    generate the decison variable matrix of offspring population
%   ParameterSet    <public>    obtain the parameter setting from user
%   ParameterFile   <public>    obtain the parameter setting from file

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        N          = 100;               % Population size
        M;                              % Number of objectives
        D;                              % Number of variables
        lower;                          % Lower bound of each decision variable
        upper;                          % Upper bound of each decision variable
        evaluation = 10000;             % Maximum number of evaluations
        operator   = @EAreal;         	% Operator function
    end
    properties(SetAccess = ?INDIVIDUAL)
        evaluated  = 0;                 % Number of evaluated individuals
    end
    properties(SetAccess = private)
        algorithm  = @NSGAII;       	% Algorithm function
        problem    = @DTLZ2;            % Problem function
        gen;                            % Current generation
        maxgen;                         % Maximum generation
        run        = 1;                 % Run No.
        runtime    = 0;                 % Runtime
        result     = {};                % Set of recorded populations
        PF;                             % True Pareto front
        parameter  = struct();      	% Parameters of functions specified by users
        mode       = 1;                 % Run mode (1.show result 2.save result 3.run outputFcn)
        outputFcn  = @GLOBAL.show;  	% Function invoked after each generation when mode = 3
    end
    methods
        %% Constructor
        function obj = GLOBAL(varargin)
        %GLOBAL - Constructor of GLOBAL class
        %
        %   GLOBAL('-Name',Value,'-Name',Value,...) returns an object with
        %   the properties specified by the inputs.
        %
        %   Example:
        %       GLOBAL('-algorithm',@NSGAII,'-problem',@DTLZ2,'-N',100,...
        %              '-M',2,'-D',10)
        
            % Initialise the parameters which can be specified by user
            obj.Set();
            proStr = {'N','M','D','algorithm','problem','operator','evaluation','run','mode','outputFcn'};
            if nargin > 0
                % The parameter setting of the environment
                IsString = find(cellfun(@ischar,varargin(1:end-1)));
                [~,Loc]  = ismember(varargin(IsString),cellfun(@(S)['-',S],proStr,'UniformOutput',false));
                for i = find(Loc)
                    obj.(varargin{IsString(i)}(2:end)) = varargin{IsString(i)+1};
                end
                % The parameter setting of the algorithm, problem and
                % operator
                MatchString = regexp(varargin(IsString),'^\-.+_parameter$');
                Loc         = cell2mat(cellfun(@(S)~isempty(S),MatchString,'UniformOutput',false));
                for i = find(Loc)
                    obj.parameter.(varargin{IsString(i)}(2:end-10)) = varargin{IsString(i)+1};
                end
            end
            % Initialise other parameters via the test problem
            obj.Set(false);
            obj.problem('init',obj,1);
            obj.Set(true);
            % Determine the run mode
            if obj.mode == 1
                obj.outputFcn = @GLOBAL.show;
            elseif obj.mode == 2
                obj.outputFcn = @GLOBAL.save;
            end
            % Add the folders of the algorithm, problem and operator to the
            % top of the search path
            addpath(fileparts(which(func2str(obj.operator))));
            addpath(fileparts(which(func2str(obj.problem))));
            addpath(fileparts(which(func2str(obj.algorithm))));
            % Register this object with INDIVIDUAL class
            INDIVIDUAL.register(obj);
        end
        %% Start running the algorithm
        function Start(obj)
        %Start - Start running the algorithm
        %
        %   obj.Start() runs the algorithm with the defined setting. This
        %   method of one GLOBAL object can only be invoked once.
        %
        %   Example:
        %       obj.Start()

            if obj.evaluated <= 0
                obj.PF = obj.problem('PF',obj,10000);
                try
                    tic;
                    obj.algorithm(obj);
                catch err
                    if strcmp(err.identifier,'GLOBAL:Termination')
                        return;
                    else
                        rethrow(err);
                    end
                end
                obj.evaluated = max(obj.evaluated,obj.evaluation);
                if isempty(obj.result)
                    obj.result = {obj.evaluated,INDIVIDUAL()};
                end
            	obj.outputFcn(obj);
            end
        end
        %% Generate an initial random population
        function Population = Initialization(obj,N)
        %Initialization - Generate an initial random population
        %
        %   P = obj.Initialization() returns an initial random population
        %   (an array of INDIVIDUAL objects) with size obj.N.
        %
        %   P = obj.Initialization(N) returns an initial random population
        %   with size N.
        %
        %   Example:
        %       P = obj.Initialization()
        
            if nargin < 2
                N = obj.N;
            end
            Population = INDIVIDUAL(obj.problem('init',obj,N));
        end
        %% Stop the running if the number of evaluations has exceeded
        function notermination = NotTermination(obj,Population)
        %NotTermination - Terminate the running if the number of
        %evaluations has exceeded
        %
        %   obj.NotTermination(P) stores the population P as the final
        %   solution set for output, and returns true if the number of
        %   evaluations has not exceeded (otherwise returns false). If the
        %   number of evaluations has exceeded, then throw an error to
        %   terminate the running forcibly.
        %
        %   This function should be invoked after each generation. The
        %   function obj.outputFcn will be invoked when obj.NotTermination
        %   has been invoked. obj.outputFcn is equal to @GLOBAL.show when
        %   obj.mode = 1, and obj.outputFcn is equal to @GLOBAL.save when
        %   obj.mode = 2.
        %
        %   The runtime of running this function and obj.outputFcn is not
        %   counted as part of the runtime of the algorithm.
        %
        %   Example:
        %       obj.NotTermination(Population)
        
            % Accumulate the runtime
            obj.runtime = obj.runtime + toc;
            % Record the last population
            index = max(1,min(min(10,size(obj.result,1)+1),ceil(10*obj.evaluated/obj.evaluation)));
            obj.result(index,:) = {obj.evaluated,Population};
            % Invoke the interrupt function
            drawnow();
            obj.outputFcn(obj);
            % Detect whether the number of evaluations has exceeded
            notermination = obj.evaluated < obj.evaluation;
            assert(notermination,'GLOBAL:Termination','Algorithm has terminated');
            tic;
        end
        %% Generate offspring population
        function Offspring = Variation(obj,Parent,nOffspring,operator,para)
        %Variation - Generate offspring population
        %
        %   O = obj.Variation(P) generates the offspring population by the
        %   default operator function obj.operator with parents P. And the
        %   number of evaluations obj.evaluated will be increased by
        %   length(O).
        %
        %   O = obj.Variation(P,N) returns only the first N offsprings at
        %   most.
        %
        %   O = obj.Variation(P,N,Fcn) generates the offspring population
        %   by the operator function Fcn.
        %
        %   O = obj.Variation(P,N,Fcn,Para) generates the offspring
        %   population by Fcn with the parameters specified by Para.
        %
        %   Example:
        %       Offspring = obj.Variation(Population)
        %       Offspring = obj.Variation(Population,1,@DE,{1,1,1,20})

            % Save the original parameters of the operator
            if nargin > 4
                field = func2str(operator);
                if isfield(obj.parameter,field)
                    oldPara = obj.parameter.(field);
                    specified                        = cellfun(@(S)~isempty(S),para);
                    obj.parameter.(field)(specified) = para(specified);
                else
                    oldPara = [];
                    obj.parameter.(field) = para;
                end
            end
            % Generate offsprings
            if nargin > 3
                Offspring = operator(obj,Parent);
            else
                Offspring = obj.operator(obj,Parent);
            end
            % Truncate the offsprings
            if nargin > 2
                obj.evaluated = obj.evaluated - max(0,length(Offspring)-nOffspring);
                Offspring     = Offspring(1:min(length(Offspring),nOffspring));
            end
            % Recover the original parameters of the operator
            if nargin > 4
                if ~isempty(oldPara)
                    obj.parameter.(field) = oldPara;
                else
                    obj.parameter = rmfield(obj.parameter,field);
                end
            end
        end
        %% Generate the decison variable matrix of offspring population
        function OffspringDec = VariationDec(obj,varargin)
        %Variation - Generate the decison variable matrix of offspring
        %population
        %
        %   Dec = obj.Variation(PDec,...) has the same function to
        %   obj.Variation(), except that the input is the decision variable
        %   matrix of the parents and the output is the decision variable
        %   matrix of the offsprings. Besides, the number of evaluations
        %   obj.evaluated will not be increased.
        %
        %   Use O = INDIVIDUAL(Dec) to obtain the offspring objects.
        %
        %   Note that it may lead to errors when the operator visits the
        %   objective values, constraint values or additional properties of
        %   the parents, since the parent objects contain only decision
        %   variables.
        %
        %   Example:
        %       Dec = obj.VariationDec(Population)
        %       Dec = obj.VariationDec(Population,1,@DE,{1,1,1,20})
        
            % Original number of evaluations
            originEvaluation = obj.evaluated;
            % Generate parent objects with only decision variables
            INDIVIDUAL.register([]);
            varargin{1} = INDIVIDUAL(varargin{1});
            % Generate offsprings
            Offspring = obj.Variation(varargin{:});
            INDIVIDUAL.register(obj);
            obj.evaluated = originEvaluation;
            % Get the decision variable matrix of offsprings
            OffspringDec = Offspring.decs;
            % Set the infeasible decision variables to boundary values
            if ~isempty(obj.lower) && ~isempty(obj.upper)
                Lower        = repmat(obj.lower,size(OffspringDec,1),1);
                Upper        = repmat(obj.upper,size(OffspringDec,1),1);
                OffspringDec = max(min(OffspringDec,Upper),Lower);
            end
        end
        %% Obtain the parameter settings from user
        function varargout = ParameterSet(obj,varargin)
        %ParameterSet - Obtain the parameter setting from user
        %
        %   [p1,p2,...] = obj.ParameterSet(v1,v2,...) returns the setting
        %   values of p1, p2, ..., where v1, v2, ... are their default
        %   values. The setting values are specified by the user with the
        %   following form:
        %   MOEA(...,'-X_parameter',{p1,p2,...},...), where X is the
        %   function name of the caller.
        %
        %   MOEA(...,'-X_parameter',{[],p2,...},...) means that parameter
        %   p1 is not specified by the user, so that it equals to its
        %   default value v1.
        %
        %   Example:
        %       [p1,p2,p3] = obj.ParameterSet(1,2,3)

            CallStack = dbstack();
            caller    = CallStack(2).name;
            varargout = varargin;
            if isfield(obj.parameter,caller)
                specified = cellfun(@(S)~isempty(S),obj.parameter.(caller));
                varargout(specified) = obj.parameter.(caller)(specified);
            end
        end
        %% Obtain the parameter settings from file
        function varargout = ParameterFile(obj,fileName,varargin)
        %ParameterFile - Obtain the parameter setting from file
        %
        %   [p1,p2,...] = obj.ParameterFile(F,v1,v2,...) returns the
        %   setting values of p1, p2, ..., where v1, v2, ... are their
        %   default values. The setting values are saved in the file F,
        %   which is a .mat file in the same folder to the caller. If the
        %   file F does not exist, then create F and save v1, v2, ... to
        %   the file.
        %
        %   Example:
        %       [p1,p2,p3] = obj.ParameterFile('Data',1,2,3)

            CallStack = dbstack('-completenames');
            fileName  = fullfile(fileparts(CallStack(2).file),[fileName,'.mat']);
            if exist(fileName,'file') == 2
                load(fileName,'Parameter');
                varargout = Parameter;
            else
                varargout = varargin;
                Parameter = varargout;
                save(fileName,'Parameter');
            end
        end
        %% Data constraint
        function set.N(obj,value)
            if ~isempty(value)
                obj.Validation(value,'int','size of population ''-N''',1);
                obj.N = value;
            end
        end
        function set.M(obj,value)
            obj.M = obj.Set('M',value,'int','number of objectives ''-M''',2);
        end
        function set.D(obj,value)
            obj.D = obj.Set('D',value,'int','number of variables ''-D''',1);
        end
        function set.lower(obj,value)
            obj.lower = obj.Set('lower',value,'');
        end
        function set.upper(obj,value)
            obj.upper = obj.Set('upper',value,'');
        end
        function set.algorithm(obj,value)
            obj.algorithm = obj.Set('algorithm',value,'function','tested algorithm ''-algorithm''');
        end
        function set.problem(obj,value)
            obj.problem = obj.Set('problem',value,'function','test problem ''-problem''');
        end
        function set.operator(obj,value)
            obj.operator = obj.Set('operator',value,'function','employed operator ''-operator''');
        end
        function set.evaluation(obj,value)
            obj.evaluation = obj.Set('evaluation',value,'int','number of evaluations ''-evaluation''',1);
        end
        function set.run(obj,value)
            obj.run = obj.Set('run',value,'int','run No. ''-run''',1);
        end
        function set.mode(obj,value)
            obj.mode = obj.Set('mode',value,'int','run mode ''-mode''',1,3);
        end
        %% Data dependence
        function value = get.gen(obj)
            value = ceil(obj.evaluated/obj.N);
        end
        function value = get.maxgen(obj)
            value = ceil(obj.evaluation/obj.N);
        end
    end
    methods(Access = private)
        %% Access the setting of public properties
        function value = Set(obj,property,value,varargin)
            persistent lock
            switch nargin
                case 1
                    lock = struct('lock',true);
                case 2
                    lock.lock = property;
                otherwise
                    if isempty(value)
                        value = obj.(property);
                    else
                        if isfield(lock,property) && (lock.(property)||lock.lock)
                            value = obj.(property);
                        else
                            obj.Validation(value,varargin{:});
                        end
                        lock.(property) = lock.lock;
                    end
            end
        end
        %% Check the validity of the specific variable
        function Validation(obj,value,Type,str,varargin)
            switch Type
                case 'folder'
                    assert(isa(value,'char') && size(value,1)==1 && size(value,2)>0,'INPUT ERROR: the %s must be a string',str);
                    assert(isdir(value),'INPUT ERROR: the folder <%s> does not exist',value);
                case 'function'
                    assert(isa(value,'function_handle'),'INPUT ERROR: the %s must be a function handle',str);
                    assert(~isempty(which(func2str(value))),'INPUT ERROR: the function <%s> does not exist',func2str(value));
                case 'logical'
                    assert(isa(value,'logical') && isscalar(value),'INPUT ERROR: the %s must be a logical scalar',str);
                case 'multichoice'
                    assert(isa(value,'char') && size(value,1)==1 && size(value,2)>0,'INPUT ERROR: the %s must be a string',str);
                    assert(ismember(value,varargin),'INPUT ERROR: the %s must be one of %s',str,strjoin(cellfun(@(S)['''',S,''''],varargin,'UniformOutput',false),' or '));
                case 'int'
                    assert(isa(value,'double') && isreal(value) && isscalar(value) && value==fix(value),'INPUT ERROR: the %s must be a integer scalar',str);
                    if ~isempty(varargin); assert(value>=varargin{1},'INPUT ERROR: the %s must be not less than %d',str,varargin{1}); end
                    if length(varargin) > 1; assert(value<=varargin{2},'INPUT ERROR: the %s must be not more than %d',str,varargin{2}); end
                    if length(varargin) > 2; assert(mod(value,varargin{3})==0,'INPUT ERROR: the %s must be a multiple of %d',str,varargin{3}); end
                case 'real'
                    assert(isa(value,'double') && isreal(value) && isscalar(value),'INPUT ERROR: the %s must be a real scalar',str);
                    if ~isempty(varargin); assert(value>=varargin{1},'INPUT ERROR: the %s must be not less than %d',str,varargin{1}); end
                    if length(varargin) > 1; assert(value<=varargin{2},'INPUT ERROR: the %s must be not more than %d',str,varargin{2}); end
            end
        end
    end
    methods(Access = private, Static)
        %% Save the result after the algorithm has been terminated
        function save(obj)
            clc; fprintf('%s on %s, %d objectives, run %d (%6.2f%%), %.2fs passed...\n',...
                         func2str(obj.algorithm),func2str(obj.problem),obj.M,obj.run,obj.evaluated/obj.evaluation*100,obj.runtime);
            if obj.evaluated >= obj.evaluation
                folder = fullfile('Data',func2str(obj.algorithm));
                [~,~]  = mkdir(folder);
                Population     = obj.result{end};
                PF             = obj.PF;
                metric.runtime = obj.runtime;
                save(fullfile(folder,sprintf('%s_%s_M%d_%d.mat',func2str(obj.algorithm),func2str(obj.problem),obj.M,obj.run)),'Population','PF','metric');
            end
        end
        %% Display the result after the algorithm has been terminated
        function show(obj)
            clc; fprintf('%s on %s, %d objectives, run %d (%6.2f%%), %.2fs passed...\n',...
                         func2str(obj.algorithm),func2str(obj.problem),obj.M,obj.run,obj.evaluated/obj.evaluation*100,obj.runtime);
            if obj.evaluated >= obj.evaluation
                % Identify the feasible and non-dominated solutions in the
                % final population
                Feasible     = find(all(obj.result{end}.cons<=0,2));
                NonDominated = NDSort(obj.result{end}(Feasible).objs,1) == 1;
                Population   = obj.result{end}(Feasible(NonDominated));
                % Calculate the metric values
                Metrics   = {@IGD};
                Score     = cellfun(@(S)GLOBAL.Metric(S,Population,obj.PF),Metrics,'UniformOutput',false);
                MetricStr = cellfun(@(S)[func2str(S),' : %.4e  '],Metrics,'UniformOutput',false);
                % Display the results
                figure('NumberTitle','off','UserData',struct(),...
                       'Name',sprintf([MetricStr{:},'Runtime : %.2fs'],Score{:},obj.runtime));
                title(sprintf('%s on %s',func2str(obj.algorithm),func2str(obj.problem)),'Interpreter','none');
                Draw(Population.objs);
                % Add new menus to the figure
                top = uimenu(gcf,'Label','Data Source');
                uimenu(top,'Label','Result (PF)',     'CallBack',{@(hObject,~,obj,P)eval('cla;Draw(P.objs);GLOBAL.cb_menu(hObject);'),obj,Population},'Checked','on');
                uimenu(top,'Label','Result (PS)',     'CallBack',{@(hObject,~,obj,P)eval('cla;Draw(P.decs);GLOBAL.cb_menu(hObject);'),obj,Population});
                uimenu(top,'Label','Result (Special)','CallBack',{@(hObject,~,obj,P)eval('obj.problem(''draw'',obj,P.decs);GLOBAL.cb_menu(hObject);'),obj,Population});
                uimenu(top,'Label','True PF',         'CallBack',{@(hObject,~,obj)eval('cla;Draw(obj.PF);GLOBAL.cb_menu(hObject);'),obj},'Separator','on');
                uimenu(top,'Label','IGD',             'CallBack',{@GLOBAL.cb_metric,obj,@IGD},'Separator','on');
                uimenu(top,'Label','Hypervolume',     'CallBack',{@GLOBAL.cb_metric,obj,@HV});
                uimenu(top,'Label','GD',              'CallBack',{@GLOBAL.cb_metric,obj,@GD});
                uimenu(top,'Label','Spacing',         'CallBack',{@GLOBAL.cb_metric,obj,@Spacing});
            end
        end
        function cb_metric(hObject,eventdata,obj,metric)
            metricName   = func2str(metric);
            MetricValues = get(gcbf,'UserData');
            % Calculate the specified metric value of each population
            if ~isfield(MetricValues,func2str(metric)) 
                tempText = text('Units','normalized','Position',[.4 .5 0],'String','Please wait ... ...'); drawnow();
                MetricValues.(metricName)(:,1) = obj.result(:,1);
                MetricValues.(metricName)(:,2) = cellfun(@(S)GLOBAL.Metric(metric,S,obj.PF),obj.result(:,2),'UniformOutput',false);
                set(gcbf,'UserData',MetricValues);
                delete(tempText);
            end
            % Display the results
            cla; Draw(cell2mat(MetricValues.(metricName)),'-k.');
            xlabel('Number of Evaluations');
            ylabel(metricName);
            GLOBAL.cb_menu(hObject);
        end
        function cb_menu(hObject)
            % Switch the selected menu
            set(get(get(hObject,'Parent'),'Children'),'Checked','off');
            set(hObject,'Checked','on');
        end
        function value = Metric(metric,Population,PF)
            % Calculate the metric value of the population
            Feasible     = find(all(Population.cons<=0,2));
            NonDominated = NDSort(Population(Feasible).objs,1) == 1;
            try
                value = metric(Population(Feasible(NonDominated)).objs,PF);
            catch
                value = NaN;
            end
        end
    end
end