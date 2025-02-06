classdef BBOB_F21 < PROBLEM
% <2009> <single> <real> <expensive/none>
% Gallagher's Gaussian 101-me peaks function

%------------------------------- Reference --------------------------------
% N. Hansen, S. Finck, R. Ros, and A. Auger. Real-parameter black-box
% optimization benchmarking 2009: Noiseless functions definitions. RR-6829,
% INRIA. 2009.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        xopt;   % Optimal decision variables
        C;      % Parameters
        y;      % Local optima
        R;      % Orthogonal (rotation) matrix
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.M = 1;
            if isempty(obj.D); obj.D = 30; end
            obj.xopt     = 4.2096874633/2*sign(rand(RandStream('mlfg6331_64','Seed',1),1,obj.D)-0.5);
            obj.lower    = zeros(1,obj.D) - 5;
            obj.upper    = zeros(1,obj.D) + 5;
            obj.encoding = ones(1,obj.D);
            alpha = 1000.^(2*(0:99)/99);
            alpha = [1000,alpha(randperm(RandStream('mlfg6331_64','Seed',2),end))];
            for i = 1 : length(alpha)
                obj.C{i} = diag(alpha(i).^(0.5*(0:obj.D-1)/(obj.D-1)))/alpha(i)^0.25;
            end
            obj.y = [rand(RandStream('mlfg6331_64','Seed',3))*8-4,rand(RandStream('mlfg6331_64','Seed',4),1,100)*10-5];
            obj.R = orth(rand(RandStream('mlfg6331_64','Seed',2),obj.D));
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,X)
            w = [10,1.1+8*(0:99)/99];
            for i = 1 : length(obj.C)
                temp(:,i) = w(i)*exp(-0.5*mean((X-obj.y(i))*obj.R'*obj.C{i}*obj.R.*(X-obj.y(i)),2));
            end
            PopObj = Tosz(10-max(temp,[],2)).^2 + fpen(X);
        end
    end
end

function X = Tosz(X)
    Xh       = zeros(size(X));
    Xh(X~=0) = log(abs(X(X~=0)));
    C1       = zeros(size(X)) + 5.5;
    C1(X>0)  = 10;
    C2       = zeros(size(X)) + 3.1;
    C2(X>0)  = 7.9;
    X = sign(X).*exp(Xh+0.049*(sin(C1.*Xh)+sin(C2.*Xh)));
end

function F = fpen(X)
    F = sum(max(0,abs(X)-5).^2,2);
end