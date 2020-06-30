classdef module_experiment_result < handle
%module_experiment_result - The class of the result in experiment module.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        Result;         % For storing the result of each algorithm on each test instance
    end
    properties(SetAccess = private)
        Algorithms;     % All the compared algorithms
        Problems;       % All the test instances
        nRuns;          % Total run times
        nPops;          % Number of saved populations
        Folder;         % Folder for saving results
    end
    methods
        %% Constructor
        function obj = module_experiment_result(Algorithms,Problems,nRuns,nPops,Folder)
            % Initialise the data
            obj.Algorithms = Algorithms;
            obj.Problems   = Problems;
            obj.nRuns      = nRuns;
            obj.nPops      = nPops;
            obj.Folder     = Folder;
            obj.Result     = cell(size(Problems,1),size(Algorithms,1),nRuns);
            % Get the actual values of N, M, D and evaluation for each test
            % instance, and put them into Problems(:,3). The four
            % properties are set to the values assigned by users when
            % invoking GLOBAL, while they are set to the actual values when
            % shown in the table or filename.
            obj.Problems(:,3) = obj.Problems(:,2);
            for i = 1 : size(obj.Problems,1)
                Global = GLOBAL('-N',obj.Problems{i,2}{1},'-M',obj.Problems{i,2}{2},'-D',obj.Problems{i,2}{3},'-evaluation',obj.Problems{i,2}{4},...
                                '-problem',[{str2func(obj.Problems{i,1})},obj.Problems{i,2}(5:end)]);
                obj.Problems{i,3}{1} = Global.N;
                obj.Problems{i,3}{2} = Global.M;
                obj.Problems{i,3}{3} = Global.D;
                obj.Problems{i,3}{4} = Global.evaluation;
                obj.Problems{i,3}{5} = Global.problem.PF(10000);
            end
        end
        %% Get the file name of a result
        function name = filename(obj,a,p,r)
            name = fullfile(obj.Folder,obj.Algorithms{a,1},sprintf('%s_%s_M%d_D%d_%d.mat',obj.Algorithms{a,1},obj.Problems{p,1},obj.Problems{p,3}{2},obj.Problems{p,3}{3},r));
        end
        %% Calculate the metric value of a result
        function data = metricValue(obj,p,a,metricFcn,metricName,showall)
            data = [];
            for i = 1 : size(obj.Result,3)
                if ~isempty(obj.Result{p,a,i})
                    if ~isfield(obj.Result{p,a,i}.metric,metricName)
                        obj.Result{p,a,i}.metric.(metricName) = arrayfun(@(index)metricFcn(str2func(metricName),obj.Result{p,a,i}.result{index,2},obj.Problems{p,3}{5}),1:size(obj.Result{p,a,i}.result,1));
                        metric = obj.Result{p,a,i}.metric;
                        save(obj.filename(a,p,i),'metric','-append');
                    end
                    try
                        if showall
                            data = [data;obj.Result{p,a,i}.metric.(metricName)];
                        else
                            data = [data,obj.Result{p,a,i}.metric.(metricName)(end)];
                        end
                    catch
                    end
                end
            end
        end
    end
end