classdef CEC2020_F9 < PROBLEM
% <single> <real>
% Composition function 2

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
            obj.O = Data{9}.o;
            obj.M = 1;
            if isempty(obj.D) || obj.D < 10
                obj.D   = 5;
                obj.Mat = Data{9}.M_5;
            elseif obj.D < 15
                obj.D   = 10;
                obj.Mat = Data{9}.M_10;
            elseif obj.D < 20
                obj.D   = 15;
                obj.Mat = Data{9}.M_15;
            else
                obj.D   = 20;
                obj.Mat = Data{9}.M_20;
            end
            obj.lower    = zeros(1,obj.D) - 100;
            obj.upper    = zeros(1,obj.D) + 100;
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            lambda = [10 1e-6 10 1];
            delta  = [10 20 30 40];
            bias   = [0 100 200 300];
            func   = {@Ackley,@Elliptic,@Griewank,@Rastrigin};
            W      = zeros(size(PopDec,1),4);
            F      = zeros(size(W));
            for i = 1 : size(W,2)
                tmp    = sum((PopDec-repmat(obj.O(i,1:size(PopDec,2)),size(PopDec,1),1)).^2,2);
                W(:,i) = 1./(sqrt(tmp)+1e-10).*exp(-tmp/2/obj.D/delta(i)^2);
                F(:,i) = func{i}((PopDec-repmat(obj.O(i,1:size(PopDec,2)),size(PopDec,1),1))*obj.Mat((i-1)*obj.D+1:i*obj.D,:)');
            end
            W = W./repmat(sum(W,2),1,size(W,2));
            PopObj = 2400 + sum(W.*(repmat(lambda,size(F,1),1).*F+repmat(bias,size(F,1),1)),2);
        end
    end
end

function F = Ackley(X)
    F = -20*exp(-0.2*sqrt(mean(X.^2,2))) - exp(mean(cos(2*pi*X),2)) + 20 + exp(1);
end

function F = Elliptic(X)
    F = sum((1e6).^(repmat(0:size(X,2)-1,size(X,1),1)/(size(X,2)-1+1e-6)).*X.^2,2);
end

function F = Griewank(X)
    X = 6*X;
    F = sum(X.^2,2)/4000 - prod(cos(X)./repmat(sqrt(1:size(X,2)),size(X,1),1),2) + 1;
end

function F = Rastrigin(X)
    X = 0.0512*X;
    F = sum(X.^2-10*cos(2*pi*X)+10,2);
end