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
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        initFcn;                % function for initializing a population
        decFcn;                 % function for repairing invalid solution
        objFcn = @(x,d)sum(x);	% objective functions
        conFcn = @(x,d)0;      	% constraint functions
    end
    methods
        %% Constructor
        function obj = UserProblem(varargin)
            isStr = find(cellfun(@ischar,varargin(1:end-1))&~cellfun(@isempty,varargin(2:end)));
            for i = isStr(ismember(varargin(isStr),{'N','maxFE','encoding','lower','upper','initFcn','decFcn','objFcn','conFcn','parameter'}))
                obj.(varargin{i}) = varargin{i+1};
            end
            if ~iscell(obj.objFcn); obj.objFcn = {obj.objFcn}; end
            if ~iscell(obj.conFcn); obj.conFcn = {obj.conFcn}; end
            obj.M = length(obj.objFcn);
            obj.D = length(obj.lower);
        end
        %% Generate initial solutions
        function Population = Initialization(obj,N)
            if nargin < 2
                N = obj.N;
            end
            if ~isempty(obj.initFcn)
                if isempty(obj.parameter)
                    PopDec = obj.initFcn(N);
                else
                    PopDec = obj.initFcn(N,obj.parameter);
                end
            elseif strcmp(obj.encoding,'binary')
                PopDec = rand(N,obj.D) < 0.5;
            elseif strcmp(obj.encoding,'permutation')
                [~,PopDec] = sort(rand(N,obj.D),2);
            else
                PopDec = unifrnd(repmat(obj.lower,N,1),repmat(obj.upper,N,1));
            end
            Population  = SOLUTION(PopDec);
            obj.optimum = max(Population.objs,[],1);
        end
        %% Repair invalid solutions
        function PopDec = CalDec(obj,PopDec)
            if ~isempty(obj.decFcn)
                for i = 1 : size(PopDec,1)
                    if isempty(obj.parameter)
                        PopDec(i,:) = obj.decFcn(PopDec(i,:));
                    else
                        PopDec(i,:) = obj.decFcn(PopDec(i,:),obj.parameter);
                    end
                end
            elseif strcmp(obj.encoding,'binary')
                PopDec = round(PopDec);
            elseif strcmp(obj.encoding,'permutation')
                [~,PopDec] = sort(PopDec,2);
                [~,PopDec] = sort(PopDec,2);
            else
                PopDec = max(min(PopDec,repmat(obj.upper,size(PopDec,1),1)),repmat(obj.lower,size(PopDec,1),1));
            end
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = zeros(size(PopDec,1),length(obj.objFcn));
            for i = 1 : size(PopDec,1)
                if isempty(obj.parameter)
                    PopObj(i,:) = cellfun(@(func)func(PopDec(i,:)),obj.objFcn);
                else
                    PopObj(i,:) = cellfun(@(func)func(PopDec(i,:),obj.parameter),obj.objFcn);
                end
            end
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,PopDec)
            PopCon = zeros(size(PopDec,1),length(obj.conFcn));
            for i = 1 : size(PopDec,1)
                if isempty(obj.parameter)
                    PopCon(i,:) = cellfun(@(func)func(PopDec(i,:)),obj.conFcn);
                else
                    PopCon(i,:) = cellfun(@(func)func(PopDec(i,:),obj.parameter),obj.conFcn);
                end
            end
        end
    end
end