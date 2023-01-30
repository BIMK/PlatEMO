classdef TP10 < PROBLEM
% <multi> <real> <large/none> <constrained> <robust>
% Test problem for robust multi-objective optimization
% delta --- 0.05 --- Maximum disturbance degree
% H     ---   50 --- Number of disturbances

%------------------------------- Reference --------------------------------
% A. Gaspar-Cunha, J. Ferreira, and G. Recio, Evolutionary robustness
% analysis for multi-objective optimization: benchmark problems, Structural
% and Multidisciplinary Optimization, 2014, 49: 771-793.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        delta;      % Maximum disturbance degree
        H;          % Number of disturbances
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            [obj.delta,obj.H] = obj.ParameterSet(0.05,50);
            obj.M = 2;
            obj.D = 3;
            obj.lower    = [0,0,1];
            obj.upper    = [10,10,3];
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(~,PopDec)
            PopObj(:,1) = PopDec(:,1).*sqrt(16+PopDec(:,3).^2) + PopDec(:,2).*sqrt(1+PopDec(:,3).^2);
            PopObj(:,2) = 20*sqrt(16+PopDec(:,3).^2)./PopDec(:,1)./PopDec(:,3);
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,PopDec)
            PopCon(:,1) = 20*sqrt(16+PopDec(:,3).^2) - 100*PopDec(:,1).*PopDec(:,3);
            PopCon(:,2) = 80*sqrt(1+PopDec(:,3).^2) - 100*PopDec(:,2).*PopDec(:,3);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = [100,100];
        end
        %% Calculate the metric value
        function score = CalMetric(obj,metName,Population)
            switch metName
                case {'Mean_IGD','Mean_HV','Worst_IGD','Worst_HV'}
                    score = feval(metName,Population,obj);
                otherwise
                    score = feval(metName,Population,obj.optimum);
            end
        end
        %% Perturb solutions multiple times
        function PopX = Perturb(obj,PopDec,N)
            if nargin < 3; N = obj.H; end
            Delta = repmat(obj.delta.*(obj.upper-obj.lower),N*size(PopDec,1),1);
            w     = UniformPoint(N,obj.D,'Latin');
            Dec   = 2*Delta.*w(reshape(repmat(1:end,size(PopDec,1),1),1,[]),:) + repmat(PopDec,N,1) - Delta;
            Dec   = obj.CalDec(Dec);
            PopX  = SOLUTION(Dec,obj.CalObj(Dec),obj.CalCon(Dec));
        end
    end
end