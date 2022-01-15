classdef UserProblem < PROBLEM
%UserProblem - User defined problem.
%
%   This problem provides a general framework, whose details can be defined
%   by the inputs of the constructor.
%
% All the acceptable properties:
%   encoding  	<string>            encoding scheme of decision variables
%   lower    	<vector>            lower bound of decision variables
%   upper      	<vector>            upper bound of decision variables
%   initFcn     <function handle>   function for initializing a population
%   decFcn      <function handle>   function for repairing invalid solution
%   objFcn     	<function handle>   objective functions
%   conFcn     	<function handle>   constraint functions
%   parameter   <any>               dataset of the problem

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = protected)
        initFcn = {};   	% Function for initializing a population
        decFcn  = {};    	% Function for repairing invalid solution
        objFcn  = {};     	% Objective functions
        conFcn  = {};     	% Constraint functions
    end
    methods
        %% Constructor
        function obj = UserProblem(varargin)
            isStr = find(cellfun(@ischar,varargin(1:end-1))&~cellfun(@isempty,varargin(2:end)));
            for i = isStr(ismember(varargin(isStr),{'N','D','maxFE','encoding','lower','upper','initFcn','decFcn','objFcn','conFcn','parameter'}))
                obj.(varargin{i}) = varargin{i+1};
            end
            useData   = ~isempty(obj.parameter);
            obj.lower = Str2Fcn(obj.lower,1,useData,'lower bound');
            obj.upper = Str2Fcn(obj.upper,1,useData,'upper bound');
            if isempty(obj.D); obj.D = length(obj.lower); end
            assert(~strcmp(obj.encoding,'real')||ismatrix(obj.lower)&&all(size(obj.lower)==[1,obj.D]),'the lower bound should be a 1*%d vector.',obj.D);
            assert(~strcmp(obj.encoding,'real')||ismatrix(obj.upper)&&all(size(obj.upper)==[1,obj.D]),'the upper bound should be a 1*%d vector.',obj.D);
            obj.parameter = Str2Fcn(obj.parameter,1,useData,'dataset');
            obj.initFcn   = Str2Fcn(obj.initFcn,2,useData,'initialization function');
            obj.decFcn    = Str2Fcn(obj.decFcn,3,useData,'repair function');
            if ~iscell(obj.objFcn); obj.objFcn = {obj.objFcn}; end
            obj.objFcn(cellfun(@isempty,obj.objFcn)) = [];
            for i = 1 : length(obj.objFcn)
                obj.objFcn{i} = Str2Fcn(obj.objFcn{i},3,useData,sprintf('objective function f%d',i));
            end
            if ~iscell(obj.conFcn); obj.conFcn = {obj.conFcn}; end
            obj.conFcn(cellfun(@isempty,obj.conFcn)) = [];
            for i = 1 : length(obj.conFcn)
                obj.conFcn{i} = Str2Fcn(obj.conFcn{i},3,useData,sprintf('constraint function g%d',i));
            end
            obj.M = length(obj.objFcn);
        end
        %% Generate initial solutions
        function Population = Initialization(obj,N)
            if nargin < 2; N = obj.N; end
            if ~isempty(obj.initFcn)
                Population = SOLUTION(CallFcn(obj.initFcn,N,obj.parameter,'initialization function',[N obj.D]));
            else
                Population = Initialization@PROBLEM(obj,N);
            end
            obj.optimum = max(Population.objs,[],1);
        end
        %% Repair invalid solutions
        function PopDec = CalDec(obj,PopDec)
            if ~isempty(obj.decFcn)
                for i = 1 : size(PopDec,1)
                    PopDec(i,:) = CallFcn(obj.decFcn,PopDec(i,:),obj.parameter,'repair function',[1 obj.D]);
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
                        PopObj(i,j) = CallFcn(obj.objFcn{j},PopDec(i,:),obj.parameter,sprintf('objective function f%d',j),[1 1]);
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
                        PopCon(i,j) = CallFcn(obj.conFcn{j},PopDec(i,:),obj.parameter,sprintf('constraint function g%d',j),[1 1]);
                    end
                end
            else
                PopCon = CalCon@PROBLEM(obj,PopDec);
            end
        end
    end
end

function f = Str2Fcn(f,type1,type2,name)
    if ischar(f)
        try
            if ~isempty(regexp(f,'^<.+>$','once'))
                switch type1
                    case 1
                        f = load(f(2:end-1));
                    otherwise
                        [folder,file] = fileparts(f(2:end-1));
                        addpath(folder);
                        f = str2func(file);
                end
            else
                switch type1
                    case 1
                        f = str2num(f);
                    case 2
                        if type2
                            f = str2func(['@(N,data)',f]);
                        else
                            f = str2func(['@(N)',f]);
                        end
                    case 3
                        if type2
                            f = str2func(['@(x,data)',f]);
                        else
                            f = str2func(['@(x)',f]);
                        end
                end
            end
        catch err
            err = addCause(err,MException('','Fail to define the %s',name));
            rethrow(err);
        end
    end
end

function output = CallFcn(func,input,parameter,name,presize)
    try
        if isempty(parameter)
            output = func(input);
        else
            output = func(input,parameter);
        end
        assert(ismatrix(output)&&all(size(output)==presize),'the size of its output should be %d*%d.',presize(1),presize(2));
    catch err
        err = addCause(err,MException('','The %s is invalid',name));
        rethrow(err);
    end
end