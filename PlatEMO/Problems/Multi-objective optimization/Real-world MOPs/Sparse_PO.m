classdef Sparse_PO < PROBLEM
% <multi> <real> <large/none> <expensive/none> <sparse/none>
% The portfolio optimization problem
% dataNo --- 1 --- Number of dataset

%------------------------------- Reference --------------------------------
% Y. Tian, C. Lu, X. Zhang, K. C. Tan, and Y. Jin, Solving large-scale
% multi-objective optimization problems with sparse optimal solutions via
% unsupervised neural networks, IEEE Transactions on Cybernetics, 2021,
% 51(6): 3115-3128.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% The datasets are the minutely closing prices of EUR/CHF taken from MT4
% No.   Name	Instruments   Length
% 1     1000	   1000         50
% 2     5000       5000         50

    properties(Access = private)
        Yield;
        Risk;
    end
    methods
    	%% Default settings of the problem
        function Setting(obj)
            % Load data
            dataNo    = obj.ParameterSet(1);
            CallStack = dbstack('-completenames');
            load(fullfile(fileparts(CallStack(1).file),'Dataset_PO.mat'),'Dataset');
            str       = {'data1000','data5000'};
            Data      = Dataset.(str{dataNo});
            obj.Yield = log(Data(:,2:end)) - log(Data(:,1:end-1));
            obj.Risk  = cov(obj.Yield');
            % Parameter setting
            obj.M = 2;
            obj.D = size(obj.Yield,1);
            obj.lower    = zeros(1,obj.D) - 1;
            obj.upper    = zeros(1,obj.D) + 1;
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopDec = PopDec./repmat(max(sum(abs(PopDec),2),1),1,size(PopDec,2));            
            PopObj = zeros(size(PopDec,1),2);
            for i = 1 : size(PopDec,1)
                PopObj(i,1) = PopDec(i,:)*obj.Risk*PopDec(i,:)';
                PopObj(i,2) = 1 - sum(PopDec(i,:)*obj.Yield);
            end
        end
        %% Display a population in the objective space
        function DrawObj(obj,Population)
            PopObj = Population.objs;
            Draw([PopObj(:,1),1-PopObj(:,2)],{'Risk','Return',[]});
        end
    end
end