classdef module_test_result < handle
%module_test_result - The class of the result in test module.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        playmaxvalue = 0;           % Maximum value of the play bar
        metric       = struct();    % Metric values
        metricline   = struct();    % Evolutionary trajectories of metric values
        result       = {};          % PF and PS in each generation
    end
    properties(SetAccess = private)
        name;                       % Name of this result
        setting;                    % String of information about setting
        Global;                     % The GLOBAL object for this result
    end
    methods
        %% Constructor
        function obj = module_test_result(name,Global)
            obj.name   = name;
            obj.Global = Global;
        end
        %% Set the properties of the result after the algorithm has terminated
        function finalset(obj)
            obj.metric.runtime = obj.Global.runtime;
            Parameter   = [repmat({' '},length(fields(obj.Global.parameter)),1),...
                           strcat({'<Parameters in '},fields(obj.Global.parameter),'>'),...
                           cellfun(@obj.getStr,struct2cell(obj.Global.parameter),'UniformOutput',false)]';
            obj.setting = ['Algorithm:  ',func2str(obj.Global.algorithm)
                           'Problem:  ',class(obj.Global.problem)
                           ' '
                           'N:  ',num2str(obj.Global.N)
                           'M:  ',num2str(obj.Global.M)
                           'D:  ',num2str(obj.Global.D)
                           'evaluation:  ',num2str(obj.Global.evaluated)
                           Parameter(:)];
        end
        function str = getStr(obj,S)
            str = cell(1,length(S));
            for i = 1 : length(S)
                if isa(S{i},'double') && isscalar(S{i})
                    str{i} = num2str(S{i});
                else
                    str{i} = '*';
                end
            end
            str = strjoin(str);
        end
    end
end