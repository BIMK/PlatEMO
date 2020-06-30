classdef MONRP < PROBLEM
% <problem> <Combinatorial MOP>
% The multi-objective next release problem
% m --- 100 --- Number of customers

%------------------------------- Reference --------------------------------
% Y. Zhang, M. Harman, and S. A. Mansouri, The multi-objective next release
% problem, Proceedings of the 9th Annual Conference on Genetic and
% Evolutionary Computation, 2007, 1129-1137.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
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
        %% Initialization
        function obj = MONRP()
            % Parameter setting
            m = obj.Global.ParameterSet(100);
            obj.Global.M = 2;
            if isempty(obj.Global.D)
                obj.Global.D = 100;
            end
            obj.Global.encoding = 'binary';
            % Randomly generate costs and values
            n    = obj.Global.D;
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
        %% A reference point for hypervolume calculation
        function P = PF(obj,N)
            P = [sum(obj.Cost),sum(obj.Value(:))];
        end
        %% Draw special figure
        function Draw(obj,PopDec)
            PopObj(:,1) = sum(repmat(obj.Cost,size(PopDec,1),1).*PopDec,2);
            PopObj(:,2) = sum(PopDec*obj.Value,2);
            cla; Draw(PopObj);
            xlabel('Cost'); ylabel('Satisfaction Score');
        end
    end
end