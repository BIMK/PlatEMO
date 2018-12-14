classdef mQAP < PROBLEM
% <problem> <Combinatorial MOP>
% The multi-objective quadratic assignment problem
% c --- 0 --- Correlation parameter

%------------------------------- Reference --------------------------------
% J. Knowles and D. Corne, Instance generators and test suites for the
% multiobjective quadratic assignment problem, Proceedings of the
% International Conference on Evolutionary Multi-Criterion Optimization,
% 2003, 295-310.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        a;  % Distance between each two locations
        b;  % Flow between each two facilities
    end
    methods
        function obj = mQAP()
            %% Initialization
            % Parameter setting
            c = obj.Global.ParameterSet(0);
            if isempty(obj.Global.M)
                obj.Global.M = 2;
            end
            if isempty(obj.Global.D)
                obj.Global.D = 10;
            end
            obj.Global.encoding = 'permutation';
            % Randomly generate distances and flows
            file = sprintf('mQAP-M%d-D%d-c%.4f.mat',obj.Global.M,obj.Global.D,c);
            file = fullfile(fileparts(mfilename('fullpath')),file);
            if exist(file,'file') == 2
                load(file,'a','b');
            else
                a = randi(100,obj.Global.D,obj.Global.D).*~eye(obj.Global.D);
                r = arrayfun(@(s)rand(obj.Global.D),1:obj.Global.M,'UniformOutput',false);
                if c >= 0.9999
                    [r{2:end}] = deal(r{1});
                elseif c >= 0
                    r(2:end) = cellfun(@(S)norminv(S.*normcdf(1,r{1},1-sqrt(c))+(1-S).*normcdf(0,r{1},1-sqrt(c)),r{1},1-sqrt(c)),r(2:end),'UniformOutput',false);
                elseif c > -0.9999
                    r(2:end) = cellfun(@(S)1-norminv(S.*normcdf(1,r{1},1-sqrt(-c))+(1-S).*normcdf(0,r{1},1-sqrt(-c)),r{1},1-sqrt(-c)),r(2:end),'UniformOutput',false);
                else
                    [r{2:end}] = deal(1-r{1});
                end
                b = cellfun(@(S)100*S.*~eye(obj.Global.D),r,'UniformOutput',false);
                save(file,'a','b');
            end
            obj.a = a;
            obj.b = b;
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = zeros(size(PopDec,1),obj.Global.M);
            [~,pi] = sort(PopDec,2);
            for i = 1 : obj.Global.M
                for j = 1 : size(PopDec,1)
                    PopObj(j,i) = sum(sum(obj.a.*obj.b{i}(pi(j,:),pi(j,:))));
                end
            end
        end
        %% A reference point for hypervolume calculation
        function P = PF(obj,N)
            P = zeros(1,obj.Global.M) + 10000*obj.Global.D.^2;
        end
    end
end