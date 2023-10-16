classdef CEC2017_F22 < PROBLEM
% <single> <real> <constrained>
% CEC'2017 constrained optimization benchmark problem

%------------------------------- Reference --------------------------------
% G. Wu, R. Mallipeddi, and P. N. Suganthan, Problem definitions and
% evaluation criteria for the CEC 2017 competition on constrained real-
% parameter optimization, National University of Defense Technology, China,
% 2016.
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
            load(fullfile(fileparts(CallStack(1).file),'CEC2017.mat'),'Data');
            obj.O = Data{12}.o;
            obj.M = 1;
            if isempty(obj.D) || obj.D < 30
                obj.D   = 10;
                obj.Mat = Data{12}.M_10;
            elseif obj.D < 50
                obj.D   = 30;
                obj.Mat = Data{12}.M_30;
            elseif obj.D < 100
                obj.D   = 50;
                obj.Mat = Data{12}.M_50;
            else
                obj.D   = 100;
                obj.Mat = Data{12}.M_100;
            end
            obj.lower    = zeros(1,obj.D) - 100;
            obj.upper    = zeros(1,obj.D) + 100;
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            Z = PopDec - repmat(obj.O(1:size(PopDec,2)),size(PopDec,1),1);
            Z = Z*obj.Mat';
            PopObj = sum(100*(Z(:,1:end-1).^2-Z(:,2:end)).^2+(Z(:,1:end-1)-1).^2,2);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,PopDec)
            Z = PopDec - repmat(obj.O(1:size(PopDec,2)),size(PopDec,1),1);
            Z = Z*obj.Mat';
            PopCon(:,1) = sum(Z.^2-10*cos(2*pi*Z)+10,2) - 100;
            PopCon(:,2) = sum(Z,2) - 2*size(Z,2);
            PopCon(:,3) = 5 - sum(Z,2);
        end
    end
end