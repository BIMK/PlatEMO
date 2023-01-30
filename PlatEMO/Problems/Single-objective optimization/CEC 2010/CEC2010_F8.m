classdef CEC2010_F8 < PROBLEM
% <single> <real> <constrained>
% CEC'2010 constrained optimization benchmark problem

%------------------------------- Reference --------------------------------
% R. Mallipeddi and P. N. Suganthan, Problem definitions and evaluation
% criteria for the CEC 2010 competition on constrained real-parameter
% optimization, Nanyang Technological University, Singapore, 2010.
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
            load(fullfile(fileparts(CallStack(1).file),'CEC2010.mat'),'Data');
            obj.O = Data{8}.O;
            obj.M = 1;
            if isempty(obj.D) || obj.D < 30
                obj.D   = 10;
                obj.Mat = Data{8}.M_10;
            else
                obj.D   = 30;
                obj.Mat = Data{8}.M_30;
            end
            if isempty(obj.D); obj.D = 10; end
            obj.D = min(obj.D,length(obj.O));
            obj.lower    = zeros(1,obj.D) - 140;
            obj.upper    = zeros(1,obj.D) + 140;
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            Z = 1 + PopDec - repmat(obj.O(1:size(PopDec,2)),size(PopDec,1),1);
            PopObj = sum(100*(Z(:,1:end-1).^2-Z(:,2:end)).^2+(Z(:,1:end-1)-1).^2,2);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,PopDec)
            Z = PopDec - repmat(obj.O(1:size(PopDec,2)),size(PopDec,1),1);
            Y = Z*obj.Mat;
            PopCon(:,1) = abs(mean((-Y.*sin(sqrt(abs(Y)))),2)) - 1e-4;
            PopCon(:,2) = 0.5 - exp(-0.1*sqrt(mean(Y.^2,2))) - 3*exp(mean(cos(0.1*Y),2)) + exp(1);
        end
    end
end