classdef CEC2020_F1 < PROBLEM
% <single> <real>
% Bent cigar function

%------------------------------- Reference --------------------------------
% C .T. Yue, K. V. Price, P. N. Suganthan, J. J. Liang, M. Z. Ali, B. Y.
% Qu, N. H. Awad, and P. P Biswas, Problem definitions and evaluation
% criteria for the CEC 2020 special session and competition on single
% objective bound constrained numerical optimization, Zhengzhou University,
% China and Nanyang Technological University, Singapore, 2019.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        O;      % Optimal decision vector
        Mat;	% Rotation matrix
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            CallStack = dbstack('-completenames');
            load(fullfile(fileparts(CallStack(1).file),'CEC2020.mat'),'Data');
            obj.O = Data{1}.o;
            obj.M = 1;
            if isempty(obj.D) || obj.D < 10
                obj.D   = 5;
                obj.Mat = Data{1}.M_5;
            elseif obj.D < 15
                obj.D   = 10;
                obj.Mat = Data{1}.M_10;
            elseif obj.D < 20
                obj.D   = 15;
                obj.Mat = Data{1}.M_15;
            else
                obj.D   = 20;
                obj.Mat = Data{1}.M_20;
            end
            obj.lower    = zeros(1,obj.D) - 100;
            obj.upper    = zeros(1,obj.D) + 100;
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            Z = PopDec - repmat(obj.O(1:size(PopDec,2)),size(PopDec,1),1);
            Y = Z*obj.Mat';
            PopObj = 100 + Y(:,1).^2 + 1e6*sum(Y(:,2:end).^2,2);
        end
    end
end