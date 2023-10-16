classdef CEC2020_F3 < PROBLEM
% <single> <real>
% Shifted and rotated Lunacek bi-Rastrigin function

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
            obj.O = Data{3}.o;
            obj.M = 1;
            if isempty(obj.D) || obj.D < 10
                obj.D   = 5;
                obj.Mat = Data{3}.M_5;
            elseif obj.D < 15
                obj.D   = 10;
                obj.Mat = Data{3}.M_10;
            elseif obj.D < 20
                obj.D   = 15;
                obj.Mat = Data{3}.M_15;
            else
                obj.D   = 20;
                obj.Mat = Data{3}.M_20;
            end
            obj.lower    = zeros(1,obj.D) - 100;
            obj.upper    = zeros(1,obj.D) + 100;
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            s   = 1 - 1/(2*sqrt(obj.D+20)-8.2);
            mu0 = 2.5;
            mu1 = -sqrt((mu0^2-1)/s);
            Y   = (PopDec-repmat(obj.O(1:size(PopDec,2)),size(PopDec,1),1))/10;
            tmp = 2*repmat(sign(obj.O(1:size(PopDec,2))),size(PopDec,1),1).*Y + mu0;
            Z   = (tmp-mu0)*obj.Mat';
            PopObj = 700 + min(sum((tmp-mu0).^2,2),obj.D+s*sum((tmp-mu1).^2,2)) + 10*(obj.D-sum(cos(2*pi*Z),2));
        end
    end
end