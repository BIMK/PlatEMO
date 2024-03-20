classdef UserProblem < PROBLEM
%UserProblem - User-defined problem.
%
%   This problem provides a general framework, whose details can be defined
%   by the inputs of the constructor.
%
% All the acceptable properties:
%   encoding  	<string>            encoding scheme of each decision variable (1.real 2.integer 3.label 4.binary 5.permutation)
%   lower    	<vector>            lower bound of each decision variable
%   upper      	<vector>            upper bound of each decision variable
%   initFcn     <function handle>   function for initializing solutions
%   evalFcn     <function handle>   function for evaluating solutions
%   decFcn      <function handle>   function for repairing invalid solutions
%   objFcn     	<function handle>   objective functions
%   conFcn     	<function handle>   constraint functions
%   gradFcn     <function handle>   function for calculating the gradients of objectives and constraints
%   data        <any>               data of the problem

%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = protected)
        initFcn    = {};        % Function for initializing solutions
        evalFcn    = {};        % Function for evaluating solutions
        decFcn     = {};    	% Function for repairing invalid solutions
        objFcn     = {};     	% Objective functions
        conFcn     = {};     	% Constraint functions
        gradFcn    = {};    	% Function for calculating the gradients of objectives and constraints
        data       = {};        % Data of the problem
    end
    methods
        %% Constructor
        function obj = UserProblem(varargin)
            isStr = find(cellfun(@ischar,varargin(1:end-1))&~cellfun(@isempty,varargin(2:end)));
            for i = isStr(ismember(varargin(isStr),{'N','M','D','maxFE','maxRuntime','encoding','lower','upper','initFcn','evalFcn','decFcn','objFcn','conFcn','gradFcn','data'}))
                obj.(varargin{i}) = varargin{i+1};
            end
            if isempty(obj.D)
                obj.encoding = Str2Fcn(obj.encoding,1,[],'encoding scheme');
                obj.D = length(obj.encoding);
            else
                obj.encoding = Str2Fcn(obj.encoding,1,[],'encoding scheme',obj.D);
            end
            obj.lower   = Str2Fcn(obj.lower,1,[],'lower bound',obj.D);
            obj.upper   = Str2Fcn(obj.upper,1,[],'upper bound',obj.D);
            obj.data    = Str2Fcn(obj.data,1,[],'dataset');
            obj.initFcn = Str2Fcn(obj.initFcn,2,~isempty(obj.data),'initialization function');
            obj.evalFcn = Str2Fcn(obj.evalFcn,3,~isempty(obj.data),'evaluation function');
            obj.decFcn  = Str2Fcn(obj.decFcn,3,~isempty(obj.data),'repair function');
            obj.objFcn  = Strs2Fcns(obj.objFcn,4,~isempty(obj.data),'objective function f');
            obj.conFcn  = Strs2Fcns(obj.conFcn,4,~isempty(obj.data),'constraint function g');
            obj.gradFcn = Str2Fcn(obj.gradFcn,3,~isempty(obj.data),'gradient function');
            Pop   = obj.Initialization(1);
            obj.M = length(Pop.objs);
        end
        %% Generate initial solutions
        function Population = Initialization(obj,N)
            if nargin < 2
                N = obj.N;
            end
            if ~isempty(obj.initFcn)
                Population = obj.Evaluation(CallFcn(obj.initFcn,N,obj.data,'initialization function',[N obj.D]));
            else
                Population = Initialization@PROBLEM(obj,N);
            end
            obj.optimum = max(Population.objs,[],1);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            if ~isempty(obj.evalFcn)
                for i = 1 : size(varargin{1},1)
                    [PopDec(i,:),PopObj(i,:),PopCon(i,:)] = CallFcn(obj.evalFcn,varargin{1}(i,:),obj.data,'evaluation function',[1 obj.D]);
                end
                Population = SOLUTION(PopDec,PopObj,PopCon,varargin{2:end});
                obj.FE     = obj.FE + length(Population);
            else
                Population = Evaluation@PROBLEM(obj,varargin{:});
            end
        end
        %% Repair invalid solutions
        function PopDec = CalDec(obj,PopDec)
            if ~isempty(obj.decFcn)
                for i = 1 : size(PopDec,1)
                    PopDec(i,:) = CallFcn(obj.decFcn,PopDec(i,:),obj.data,'repair function',[1 obj.D]);
                end
            else
                PopDec = CalDec@PROBLEM(obj,PopDec);
            end
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            if ~isempty(obj.objFcn)
                PopObj = zeros(size(PopDec,1),length(obj.objFcn));
                for i = 1 : size(PopDec,1)
                    for j = 1 : length(obj.objFcn)
                        PopObj(i,j) = CallFcn(obj.objFcn{j},PopDec(i,:),obj.data,sprintf('objective function f%d',j),[1 1]);
                    end
                end
            else
                PopObj = CalObj@PROBLEM(obj,PopDec);
            end
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,PopDec)
            if ~isempty(obj.conFcn)
                PopCon = zeros(size(PopDec,1),length(obj.conFcn));
                for i = 1 : size(PopDec,1)
                    for j = 1 : length(obj.conFcn)
                        PopCon(i,j) = CallFcn(obj.conFcn{j},PopDec(i,:),obj.data,sprintf('constraint function g%d',j),[1 1]);
                    end
                end
            else
                PopCon = CalCon@PROBLEM(obj,PopDec);
            end
        end
        %% Calculate the gradients of objectives and constraints
        function [ObjGrad,ConGrad] = CalGrad(obj,Dec)
            if ~isempty(obj.gradFcn)
            	[ObjGrad,ConGrad] = CallFcn(obj.gradFcn,Dec,obj.data,'gradient function');
            else
                [ObjGrad,ConGrad] = CalGrad@PROBLEM(obj,Dec);
            end
        end
    end
end

function var = Str2Fcn(var,type,useData,name,D)
% Convert a string into a variable or function and check its validity

    if ischar(var)
        try
            if ~isempty(regexp(var,'^<.+>$','once'))
                switch type
                    case 1      % For lower, upper, data
                        var = load(var(2:end-1));
                    otherwise   % For initFcn, evalFcn, decFcn, objFcn, conFcn, objGradFcn, objConFcn
                        [folder,file,ext] = fileparts(var(2:end-1));
                        if type ~= 4 || strcmp(ext,'.m')
                            addpath(folder);
                            var = str2func(file);
                        else
                            var = load(var(2:end-1));
                        end
                end
            else
                switch type
                    case 1      % For lower, upper, data
                        var = str2num(var);
                    case 2      % For initFcn
                        if useData
                            var = str2func(['@(N,data)',var]);
                        else
                            var = str2func(['@(N)',var]);
                        end
                    case {3,4}	% For evalFcn, decFcn, objFcn, conFcn, objGradFcn, objConFcn
                        if useData
                            var = str2func(['@(x,data)',var]);
                        else
                            var = str2func(['@(x)',var]);
                        end
                end
            end
        catch err
            err = addCause(err,MException('','Fail to define the %s',name));
            rethrow(err);
        end
    end
    if type == 1 && nargin > 4      % For lower, upper, data
        if isscalar(var)
            var = repmat(var,1,D);
        else
            assert(ismatrix(var)&&all(size(var)==[1,D]),'the %s should be a scalar or a 1*%d vector, while its current size is %d*%d.',name,D,size(var,1),size(var,2));
        end
    end
    if type == 4 && isnumeric(var)   % For objFcn, conFcn
        try
            fprintf('Fit the %s...\n',name);
            Model = fitrgp(var(:,1:end-1),var(:,end),'OptimizeHyperparameters','all','HyperparameterOptimizationOptions',struct('ShowPlots',false,'Verbose',0));
            var   = @(x)predict(Model,x);
        catch err
            err = addCause(err,MException('','Fail to fit the %s',name));
            rethrow(err);
        end
    end
end

function Var = Strs2Fcns(Var,type,useData,name)
% Convert multiple strings into functions

    if ~iscell(Var)
        Var = {Var};
    end
    Var(cellfun(@isempty,Var)) = [];
    for i = 1 : length(Var)
        Var{i} = Str2Fcn(Var{i},type,useData,[name,num2str(i)]);
    end
end

function varargout = CallFcn(func,input,data,name,varargin)
% Call a function and check the validity of its output

    try
        if isempty(data)
            [varargout{1:nargout}] = func(input);
        else
            [varargout{1:nargout}] = func(input,data);
        end
        for i = 1 : min(length(varargout),length(varargin))
            assert(ismatrix(varargout{i})&&all(size(varargout{i})==varargin{i}),'the size of its output #%d should be %d*%d, while its current size is %d*%d.',i,varargin{i}(1),varargin{i}(2),size(varargout{i},1),size(varargout{i},2));
        end
    catch err
        err = addCause(err,MException('','The %s is invalid',name));
        rethrow(err);
    end
end