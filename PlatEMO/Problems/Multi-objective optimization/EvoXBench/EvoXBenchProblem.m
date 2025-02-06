classdef EvoXBenchProblem < PROBLEM
%EvoXBenchProblem - The base problems of Multiobjective Neural Architecture Search.
%
%   EvoXBenchProblem is a subclass of PROBLEM. In EvoXBench, offspring 
%   generated by algorithms are always evaluated on the training set.When 
%   saving the results, it is necessary to re-evaluate on the test set, 
%   which can also obtain accurate objective values.
%   Before running the EvoXBench test suite, you need to install evoxbench, 
%   the installation tutorial address is: https://github.com/EMI-Group/evoxbench.git
%
% EvoXBenchProblem properties:
%   conn            <class>     the TCP/IP object
%   true_eval       <scalar>    true is validation pthase or test phase, false is training phase
%
% EvoXBenchProblem methods:
%   delete          <public>    release resources in Python
%   InitRemoteObj   <public>    initializes the TCP/IP object
%   Setting         <public>    default settings of the problem
%   CalDec          <public>    repair invalid solutions
%   CalObj          <public>    calculate the objective values of solutions
%   Sample          <public>    sample the network architectures
%   GetOptimum      <public>    generate the optimal objective values of the problem
%   GetPF          	<public>    generate the image of the Pareto front
%   CalMetric       <public>    calculate the metric value of a population
%   DrawObj         <public>    display solutions in the objective space
%   SetValid        <public>    set to validation phase
%   SetTest         <public>    set to test phase
%
%   Example:
%
%       platemo()
%
%   displays the GUI of PlatEMO.When users need to save the population on 
%   EvoXBench, the objective value of the population needs to be 
%   recalculated on the test set. Since the saved results are evaluated on 
%   the training set, the decision variables are accurate, but objective 
%   values are less accurate. 
%
%       pro = @C10MOP1;
%       problem = pro();
%       [decs, objs, cons] = platemo('algorithm',@NSGAII,'problem',pro,'N',100,'maxFE',10000);
%       population = SOLUTION(decs, objs, cons);
%       igd = problem.CalMetric('IGD', population);
%   runs NSGAII with a population size of 100 on C10MOP1 for 10000
%   evaluations and calculate the IGD value. If the users want to save the 
%   population, repair invalid solutions and recalculate the objective 
%   value.

%------------------------------- Reference --------------------------------
% Z. Lu, R. Cheng, Y. Jin, K. C. Tan, and K. Deb. Neural architecture
% search as multiobjective optimization benchmarks: Problem formulation and
% performance assessment. IEEE Transactions on Evolutionary Computation,
% 2024, 28(2): 323-337.
%--------------------------------------------------------------------------

    properties (Access = private)
        conn;
        true_eval;
    end
    methods
        %% Release resources in Python
        function delete(obj)
            if ~isempty(obj.conn)
                request.operation = 'delete';
                try
                    obj.conn.writeline(jsonencode(request));
                    obj.conn.flush('output');
                catch
                end
            end
        end
        %% Initializes the TCP/IP object
        function InitRemoteObj(obj, config)
            try
                request.config = config;
                request.operation = 'create';
                obj.conn = tcpclient('127.0.0.1', 9876, 'timeout', 30);
                obj.conn.writeline(jsonencode(request));
                obj.conn.flush('output');
                response = jsondecode(obj.conn.readline());
                assert(strcmp(response.status,'ok'),'Failed to connect to the remote server.');
            catch err
                error('EvoXBench:createFailed','%s\nPlease visit https://github.com/EMI-Group/evoxbench for details.',err.message);
            end
        end
        %% Default settings of the problem
        function Setting(obj, config)
            obj.InitRemoteObj(config);
            query.operation = 'settings';
            obj.conn.writeline(jsonencode(query));
            obj.conn.flush('output');
            response = jsondecode(obj.conn.readline());
            if strcmp(response.status, 'ok')
                basicinfo = response.result;
            else
                exception = MException('TestFunction:queryFailed','Failed to query');
                throw(exception);
            end
            obj.M = basicinfo.n_obj;
            if isempty(obj.D); obj.D = basicinfo.n_var; end
            obj.lower    = reshape(basicinfo.lb, 1, []);
            obj.upper    = reshape(basicinfo.ub, 1, []);
            obj.encoding = ones(1,obj.D);
            obj.true_eval = false;
        end
        %% Repair invalid solutions
        function PopDec = CalDec(obj, PopDec)
            % clip and round, because evoxbench uses discrete encoding
            PopDec = min(PopDec, repmat(obj.upper, size(PopDec,1), 1));
            PopDec = max(PopDec, repmat(obj.lower, size(PopDec,1), 1));
            PopDec = round(PopDec);

            % deduplicate
            TotalLen = size(PopDec, 1);
            [~, ia] = unique(PopDec, 'rows');
            NonDupLen = size(ia, 1);
            DupLen = TotalLen - NonDupLen;
            DupIndex = true(TotalLen, 1);
            DupIndex(ia) = false;
            PopDec(DupIndex, :) = obj.Sample(DupLen);
        end
        %% calculate the objective values
        function PopObj = CalObj(obj, x)
            query.operation = 'query';
            query.encoding = x;
            query.true_eval = obj.true_eval;
            obj.conn.writeline(jsonencode(query));
            obj.conn.flush('output');
            response = jsondecode(obj.conn.readline());
            if strcmp(response.status, 'ok')
                PopObj = response.result;
            else
                exception = MException('EvoXBench:queryFailed','Failed to query');
                throw(exception);
            end
        end
        %% Sample the network architectures
        function R = Sample(obj, n_samples)
            query.operation = 'sample';
            query.n_samples = n_samples;
            obj.conn.writeline(jsonencode(query));
            obj.conn.flush('output');
            response = jsondecode(obj.conn.readline());
            if strcmp(response.status, 'ok')
                R = response.result;
            else
                exception = MException('EvoXBench:queryFailed','Failed to query');
                throw(exception);
            end
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj, N)
            query.operation = 'pareto_front';
            query.n_samples = N;
            obj.conn.writeline(jsonencode(query));
            obj.conn.flush('output');
            response = jsondecode(obj.conn.readline());
            if strcmp(response.status, 'ok')
                R = response.result;
            else
                exception = MException('EvoXBench:queryFailed','Failed to query');
                throw(exception);
            end
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            % -1 will get all points
            R = obj.GetOptimum(-1);
        end
        function score = CalMetric(obj,metName,Population)
        %CalMetric - calculate the metric value of a population
        %
        %   value = obj.CalMetric(Met,P) returns the metric value of a
        %   population P, where Met is a string denoting the name of a
        %   metric function. If Met is 'HV', 'IGD' or 'normalized_hv', the 
        %   metric value is calculate through the Python side. Before the 
        %   metric value is calculated, the test phase is set and the 
        %   objective value is recalculated. Also, invalid solutions are 
        %   removed.
        %   
        %   Example:
        %       value = Problem.CalMetric('HV',Population);
        
            obj.true_eval = true;
            PopDec = obj.CalDec(Population.decs());
            PopObj = obj.CalObj(PopDec);
            PopCon = obj.CalCon(PopDec);
            invalidSolutions = any(PopObj == inf, 2);
            validSolutions = ~invalidSolutions;
            PopDec = PopDec(validSolutions, :);
            indicator = lower(metName);
            if strcmp(indicator, 'hv') || strcmp(indicator, 'igd') || strcmp(indicator, 'normalized_hv')
                query.operation = 'calc_perf_indicator';
                query.inputs = PopDec;
                query.indicator = indicator;
                obj.conn.writeline(jsonencode(query));
                obj.conn.flush('output');
                response = jsondecode(obj.conn.readline());
                if strcmp(response.status, 'ok')
                    score = response.result;
                else
                    exception = MException('EvoXBench:queryFailed','Failed to query');
                    throw(exception);
                end
            else
                Population = SOLUTION(PopDec,PopObj,PopCon);
                score = feval(metName,Population,obj.optimum);
            end
        end
        function DrawObj(obj,Population)
            %DrawObj - Display a population in the objective space.
            %
            %   To display the results of test phase in real time in the 
            %   GUI, set to validation phase. Then recalculate the 
            %   objective value and create a new SOLUTION object. 
            %   
            %	This function is usually called by the GUI.
            obj.true_eval = true;
            PopDec = obj.CalDec(Population.decs());
            PopObj = obj.CalObj(PopDec);
            PopCon = obj.CalCon(PopDec);
            Population = SOLUTION(PopDec,PopObj,PopCon);

            ax = Draw(Population.objs,{'\it f\rm_1','\it f\rm_2','\it f\rm_3'});
            if ~isempty(obj.PF)
                if ~iscell(obj.PF)
                    if obj.M == 2
                        plot(ax,obj.PF(:,1),obj.PF(:,2),'-k','LineWidth',1);
                    elseif obj.M == 3
                        plot3(ax,obj.PF(:,1),obj.PF(:,2),obj.PF(:,3),'-k','LineWidth',1);
                    end
                else
                    if obj.M == 2
                        surf(ax,obj.PF{1},obj.PF{2},obj.PF{3},'EdgeColor','none','FaceColor',[.85 .85 .85]);
                    elseif obj.M == 3
                        surf(ax,obj.PF{1},obj.PF{2},obj.PF{3},'EdgeColor',[.8 .8 .8],'FaceColor','none');
                    end
                    set(ax,'Children',ax.Children(flip(1:end)));
                end
            elseif size(obj.optimum,1) > 1 && obj.M < 4
                if obj.M == 2
                    plot(ax,obj.optimum(:,1),obj.optimum(:,2),'.k');
                elseif obj.M == 3
                    plot3(ax,obj.optimum(:,1),obj.optimum(:,2),obj.optimum(:,3),'.k');
                end
            end
            obj.true_eval = false;
        end
        %% Set to validation phase
        function SetValid(obj)
            obj.true_eval = true;
        end
        %% Set to test phase
        function SetTest(obj)
            obj.true_eval = true;
        end
    end
end