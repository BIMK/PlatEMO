classdef mQAP < PROBLEM
% <multi/many> <permutation> <large/none>
% The multi-objective quadratic assignment problem
% c --- 0 --- Correlation parameter

%------------------------------- Reference --------------------------------
% J. Knowles and D. Corne, Instance generators and test suites for the
% multiobjective quadratic assignment problem, Proceedings of the
% International Conference on Evolutionary Multi-Criterion Optimization,
% 2003, 295-310.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        a;	% Distance between each two locations
        b;	% Flow between each two facilities
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            % Parameter setting
            c = obj.ParameterSet(0);
            if isempty(obj.M); obj.M = 2; end
            if isempty(obj.D); obj.D = 10; end
            obj.encoding = 5 + zeros(1,obj.D);
            % Randomly generate distances and flows
            file = sprintf('mQAP-M%d-D%d-c%.4f.mat',obj.M,obj.D,c);
            file = fullfile(fileparts(mfilename('fullpath')),file);
            if exist(file,'file') == 2
                load(file,'a','b');
            else
                a = randi(100,obj.D,obj.D).*~eye(obj.D);
                r = arrayfun(@(s)rand(obj.D),1:obj.M,'UniformOutput',false);
                if c >= 0.9999
                    [r{2:end}] = deal(r{1});
                elseif c >= 0
                    r(2:end) = cellfun(@(S)norminv(S.*normcdf(1,r{1},1-sqrt(c))+(1-S).*normcdf(0,r{1},1-sqrt(c)),r{1},1-sqrt(c)),r(2:end),'UniformOutput',false);
                elseif c > -0.9999
                    r(2:end) = cellfun(@(S)1-norminv(S.*normcdf(1,r{1},1-sqrt(-c))+(1-S).*normcdf(0,r{1},1-sqrt(-c)),r{1},1-sqrt(-c)),r(2:end),'UniformOutput',false);
                else
                    [r{2:end}] = deal(1-r{1});
                end
                b = cellfun(@(S)100*S.*~eye(obj.D),r,'UniformOutput',false);
                save(file,'a','b');
            end
            obj.a = a;
            obj.b = b;
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = zeros(size(PopDec,1),obj.M);
            [~,pi] = sort(PopDec,2);
            for i = 1 : obj.M
                for j = 1 : size(PopDec,1)
                    PopObj(j,i) = sum(sum(obj.a.*obj.b{i}(pi(j,:),pi(j,:))));
                end
            end
        end
        %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = zeros(1,obj.M) + 10000*obj.D.^2;
        end
    end
end