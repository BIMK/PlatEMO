classdef ZXH_CF9 < PROBLEM
% <multi/many> <real> <large/none> <constrained>
% Constrained benchmark MOP proposed by Zhou, Xiang, and He

%------------------------------- Reference --------------------------------
% Y. Zhou, Y. Xiang, and X. He, Constrained multiobjective optimization:
% Test problem construction and performance evaluations, IEEE Transactions
% on Evolutionary Computation, 2021, 25(1): 172-186.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        k;  % Number of constrained variables
    end
	methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 3; end
            if isempty(obj.D); obj.D = obj.M+10;  end
            obj.lower    = zeros(1,obj.D) + 1e-10;
            obj.upper    = ones(1,obj.D)  - 1e-10;
            obj.encoding = ones(1,obj.D);
            if obj.M <= 3
                obj.k = obj.M - 1;
            elseif obj.M > 3 && obj.M <= 8 
                obj.k = floor(obj.M/2); 
            else
                obj.k = 3; 
            end
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            PopDec = varargin{1};
            OptX   = 0.2;
            [N,D]  = size(PopDec);
            M      = obj.M;
            % Step 1: Compute cumsum 
            Sx = cumsum(PopDec(:,1:M).^2,2,'reverse'); 
            % Step 2: Compute theta
            THETA = 2/pi*atan(sqrt(Sx(:,2:end))./PopDec(:,1:M-1));
            % Step 3: Calculate Ackley function
            h = 20 - 20 * exp(-0.2 * sqrt(sum((PopDec(:,M+1:end)-OptX).^2,2)/(D-M))) + exp (1) - exp(sum(cos(2 * pi .*(PopDec(:,M+1:end)-OptX)),2)/(D-M));
            % Step 4: Compute T_
            T = (1 - Sx(:,1)).^2 + h;
            % Step 5: Objectives (mixed)
            A      = 2; % number of segments
            G      = 1-[ones(N,1) cumprod(sin(pi/2*THETA),2)] .* [cos(pi/2*THETA) ones(N,1)];
            G(:,1) = THETA(:,1) - cos(2*pi*A*THETA(:,1)+pi/2)/2/A/pi;
            PopObj = G .* repmat((1+T),1,M);
            % Step 6: Constraints
            PopCon(:,1) = Sx(:,1) + h - 1; 
            for i = 1 : obj.k
                PopCon(:,i+1) = min(THETA(:,i)-1/4,3/4-THETA(:,i));
            end
            Population = SOLUTION(varargin{1},PopObj,PopCon,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,obj.M);
            R = 1 - R./repmat(sqrt(sum(R.^2,2)),1,obj.M);
            V = asin(sqrt(sum((1-R(:,2:end)).^2,2)));
            R(:,1) = (2/pi).*(V+1/8*sin(8*V));
            T = zeros(size(R));
            for i = obj.M-1 : -1 : 2
                T(:,i) = atan((1-R(:,i+1))./(1-R(:,i))./cos(T(:,i+1)));
            end
            T(:,1) = asin((1-R(:,2))./cos(T(:,2)));
            THETA  = T(:,1:obj.k)*2/pi;
            Valid  = all(THETA<=1/4|THETA>=3/4,2);
            R      = R(Valid,:);
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                R = obj.GetOptimum(100);
            elseif obj.M == 3
                a = linspace(0,pi/2,30)';
                R = {1-sin(a)*cos(a'),1-sin(a)*sin(a'),1-cos(a)*ones(size(a'))};
                V = asin(sqrt((1-R{2}).^2+(1-R{3}).^2));
                R{1} = (2/pi).*(V+1/8*sin(8*V));
                T2 = atan((1-R{3})./(1-R{2}));
                T1 = asin((1-R{2})./cos(T2));
                THETA = cat(3,T1,T2)*2/pi;
                Valid = all(THETA<=1/4|THETA>=3/4,3);
                R{1}(~Valid) = nan;
            else
                R = [];
            end
        end
    end
end