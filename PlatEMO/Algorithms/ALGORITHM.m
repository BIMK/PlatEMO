classdef ALGORITHM < handle & matlab.mixin.Heterogeneous
%ALGORITHM - The superclass of algorithms.
%
%   This is the superclass of algorithms. An object of ALGORITHM stores the
%   settings of the algorithm and the data generated in current execution.
%
% ALGORITHM properties:
%   parameter       <any>       parameters of the algorithm
%   save            <scalar>    number of populations saved in an execution
%   outputFcn       <function>	function called after each generation
%   pro             <class>     problem solved in current execution
%   result          <cell>      populations saved in current execution
%   metric          <struct>    metric values of current populations
%   starttime       <scalar>	Used for runtime recording
%
% ALGORITHM methods:
%   ALGORITHM       <protected> the constructor setting all the properties specified by user
%   Solve           <public>    use the algorithm to solve a problem
%   main            <public>	the main function of the algorithm
%   NotTerminated 	<protected>	the function called after each generation of the execution
%   ParameterSet	<protected>	obtain the parameters of the algorithm

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = protected)
        parameter = {};                 % Parameters of the algorithm
        save      = -10;            	% Number of populations saved in an execution
        outputFcn = @DefaultOutput;     % Function called after each generation
        pro;                            % Problem solved in current execution
        result;                         % Populations saved in current execution
        metric;                         % Metric values of current populations
        starttime;                      % Used for runtime recording
    end
    methods(Access = protected)
        function obj = ALGORITHM(varargin)
        %ALGORITHM - The constructor of ALGORITHM.
        %
        %   Alg = algName('Name',Value,'Name',Value,...) generates an
        %   object with the properties specified by the inputs. algName is
        %   a subclass of ALGORITHM, while ALGORITHM cannot be instantiated
        %   directly.
        %
        %   Example:
        %        Algorithm = MOEAD('parameter',4,'save',1)

            isStr = find(cellfun(@ischar,varargin(1:end-1))&~cellfun(@isempty,varargin(2:end)));
            for i = isStr(ismember(varargin(isStr),{'parameter','save','outputFcn'}))
                obj.(varargin{i}) = varargin{i+1};
            end
        end
    end
    methods(Sealed)
        function Solve(obj,Problem)
        %Solve - Use the algorithm to solve a problem.
        %
        %   obj.Solve(Pro) uses the algorithm to solve a problem, where Pro
        %   is a PROBLEM object.
        %
        %   In terms of the default obj.outputFcn, the result will be
        %   displayed when obj.save = 0 and saved when obj.save > 0.
        %
        %   Example:
        %       Algorithm.Solve(Problem)
            
            try
                obj.result = {};
                obj.metric = struct('runtime',0);
                obj.pro    = Problem;
                obj.pro.FE = 0;
                addpath(fileparts(which(class(obj))));
                addpath(fileparts(which(class(obj.pro))));
                obj.starttime = tic;
                obj.main(obj.pro);
            catch err
                if ~strcmp(err.identifier,'PlatEMO:Termination')
                    rethrow(err);
                end
            end
        end
    end
    methods
        function main(obj,Problem)
        %main - The main function of the algorithm.
        %
        %   This function is expected to be implemented in each subclass of
        %   ALGORITHM, which is usually called by ALGORITHM.Solve.
        end
    end
    methods(Access = protected, Sealed)
        function nofinish = NotTerminated(obj,Population)
        %NotTerminated - The function called after each generation of the
        %execution.
        %
        %   obj.NotTerminated(P) stores the population P as the result of
        %   the current execution, and returns true if the algorithm should
        %   be terminated, i.e., the number of function evaluations or
        %   runtime exceeds.
        %
        %   obj.outputFcn is called here, whose runtime will not be counted
        %   in the runtime of current execution.
        %
        %   Example:
        %       while Algorithm.NotTerminated(Population)
        %           ... ...
        %       end
        
            obj.metric.runtime = obj.metric.runtime + toc(obj.starttime);
            if obj.pro.maxRuntime < inf
                obj.pro.maxFE = obj.pro.FE*obj.pro.maxRuntime/obj.metric.runtime;
            end
            num   = max(1,abs(obj.save));
            index = max(1,min(min(num,size(obj.result,1)+1),ceil(num*obj.pro.FE/obj.pro.maxFE)));
            obj.result(index,:) = {obj.pro.FE,Population};
            drawnow('limitrate');
            obj.outputFcn(obj,obj.pro);
            nofinish = obj.pro.FE < obj.pro.maxFE;
            assert(nofinish,'PlatEMO:Termination','');
            obj.starttime = tic;
        end
        function varargout = ParameterSet(obj,varargin)
        %ParameterSet - Obtain the parameters of the algorithm.
        %
        %   [p1,p2,...] = obj.ParameterSet(v1,v2,...) sets the values of
        %   parameters p1, p2, ..., where each parameter is set to the
        %   value given in obj.parameter if obj.parameter is specified, and
        %   set to the value given in v1, v2, ... otherwise.
        %
        %   Example:
        %       [p1,p2,p3] = Algorithm.ParameterSet(1,2,3)

            varargout = varargin;
            specified = ~cellfun(@isempty,obj.parameter);
            varargout(specified) = obj.parameter(specified);
        end
    end
    methods(Sealed)
        function Scores = CalMetric(obj,metName)
        % Calculate metric values
        
            if ~isfield(obj.metric,metName)
                obj.metric.(metName) = [cell2mat(obj.result(:,1)),cellfun(@(S)obj.pro.CalMetric(metName,S),obj.result(:,2))];
            end
            Scores = obj.metric.(metName);
        end
    end
end

function DefaultOutput(Algorithm,Problem)
% The default output function of ALGORITHM

    clc; fprintf('%s on %d-objective %d-variable %s (%6.2f%%), %.2fs passed...\n',class(Algorithm),Problem.M,Problem.D,class(Problem),Problem.FE/Problem.maxFE*100,Algorithm.metric.runtime);
    if Problem.FE >= Problem.maxFE
        if Algorithm.save < 0
            if Problem.M > 1
                Population = Algorithm.result{end};
                if length(Population) >= size(Problem.optimum,1); name = 'HV'; else; name = 'IGD'; end
                value = Algorithm.CalMetric(name);
                figure('NumberTitle','off','Name',sprintf('%s : %.4e  Runtime : %.2fs',name,value(end),Algorithm.CalMetric('runtime')));
                title(sprintf('%s on %s',class(Algorithm),class(Problem)),'Interpreter','none');
                top = uimenu(gcf,'Label','Data source');
                g   = uimenu(top,'Label','Population (obj.)','CallBack',{@(h,~,Pro,P)eval('Draw(gca);Pro.DrawObj(P);cb_menu(h);'),Problem,Population});
                uimenu(top,'Label','Population (dec.)','CallBack',{@(h,~,Pro,P)eval('Draw(gca);Pro.DrawDec(P);cb_menu(h);'),Problem,Population});
                uimenu(top,'Label','True Pareto front','CallBack',{@(h,~,P)eval('Draw(gca);Draw(P,{''\it f\rm_1'',''\it f\rm_2'',''\it f\rm_3''});cb_menu(h);'),Problem.optimum});
                cellfun(@(s)uimenu(top,'Label',s,'CallBack',{@(h,~,A)eval('Draw(gca);Draw(A.CalMetric(h.Label),''-k.'',''LineWidth'',1.5,''MarkerSize'',10,{''Number of function evaluations'',strrep(h.Label,''_'','' ''),[]});cb_menu(h);'),Algorithm}),{'IGD','HV','GD','Feasible_rate'});
                set(top.Children(4),'Separator','on');
                g.Callback{1}(g,[],Problem,Population);
            else
                best = Algorithm.CalMetric('Min_value');
                if isempty(best); best = nan; end
                figure('NumberTitle','off','Name',sprintf('Min value : %.4e  Runtime : %.2fs',best(end),Algorithm.CalMetric('runtime')));
                title(sprintf('%s on %s',class(Algorithm),class(Problem)),'Interpreter','none');
                top = uimenu(gcf,'Label','Data source');
                uimenu(top,'Label','Population (dec.)','CallBack',{@(h,~,Pro,P)eval('Draw(gca);Pro.DrawDec(P);cb_menu(h);'),Problem,Algorithm.result{end}});
                cellfun(@(s)uimenu(top,'Label',s,'CallBack',{@(h,~,A)eval('Draw(gca);Draw(A.CalMetric(h.Label),''-k.'',''LineWidth'',1.5,''MarkerSize'',10,{''Number of function evaluations'',strrep(h.Label,''_'','' ''),[]});cb_menu(h);'),Algorithm}),{'Min_value','Feasible_rate'});
                set(top.Children(2),'Separator','on');
                top.Children(2).Callback{1}(top.Children(2),[],Algorithm);
            end
        elseif Algorithm.save > 0
            folder = fullfile('Data',class(Algorithm));
            [~,~]  = mkdir(folder);
            file   = fullfile(folder,sprintf('%s_%s_M%d_D%d_',class(Algorithm),class(Problem),Problem.M,Problem.D));
            runNo  = 1;
            while exist([file,num2str(runNo),'.mat'],'file') == 2
                runNo = runNo + 1;
            end
            result = Algorithm.result;
            metric = Algorithm.metric;
            save([file,num2str(runNo),'.mat'],'result','metric');
        end
    end
end

function cb_menu(h)
% Switch between the selected menu

    set(get(get(h,'Parent'),'Children'),'Checked','off');
    set(h,'Checked','on');
end