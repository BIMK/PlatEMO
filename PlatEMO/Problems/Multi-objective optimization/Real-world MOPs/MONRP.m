classdef MONRP < PROBLEM
% <multi> <binary> <large/none>
% The multi-objective next release problem
% m --- 100 --- Number of customers

%------------------------------- Reference --------------------------------
% Y. Zhang, M. Harman, and S. A. Mansouri, The multi-objective next release
% problem, Proceedings of the Annual Conference on Genetic and Evolutionary
% Computation, 2007, 1129-1137.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        Cost;   % Cost of each requirement
        Value;  % Value of each customer on each requirement
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            % Parameter setting
            m = obj.ParameterSet(100);
            obj.M = 2;
            if isempty(obj.D); obj.D = 100; end
            obj.encoding = 4 + zeros(1,obj.D);
            % Randomly generate costs and values
            n    = obj.D;
            file = sprintf('MONRP-n%d-m%d.mat',n,m);
            file = fullfile(fileparts(mfilename('fullpath')),file);
            if exist(file,'file') == 2
                load(file,'Cost','Value');
            else
                Cost  = randi(9,1,n);
                Value = randi([0 5],n,m);
                save(file,'Cost','Value');
            end
            obj.Cost  = Cost;
            obj.Value = Value;
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj(:,1) = sum(repmat(obj.Cost,size(PopDec,1),1).*PopDec,2);
            PopObj(:,2) = sum(obj.Value(:)) - sum(PopDec*obj.Value,2);
        end
        %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [sum(obj.Cost),sum(obj.Value(:))];
        end
        %% Display a population in the objective space
        function DrawObj(obj,Population)
            PopObj(:,1) = sum(repmat(obj.Cost,length(Population),1).*Population.decs,2);
            PopObj(:,2) = sum(Population.decs*obj.Value,2);
            Draw(PopObj,{'Cost','Satisfaction score',[]});
        end
    end
end