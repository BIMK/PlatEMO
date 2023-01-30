classdef ZXH_CF4 < PROBLEM
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

    methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 3; end
            if isempty(obj.D); obj.D = obj.M+10;  end
            obj.lower    = zeros(1,obj.D) + 1e-10;
            obj.upper    = ones(1,obj.D)  - 1e-10;
            obj.encoding = ones(1,obj.D); 
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
            % Step 3: Calculate Griewank function
            h = 5*(sum((PopDec(:,M+1:end)-OptX).^2,2)-prod(cos(10*pi*(PopDec(:,M+1:end)-OptX)./repmat(sqrt(1:(D-M)),N,1)),2)+1);
            % Step 4: Compute T_
            T = (1 - Sx(:,1)).^2 + h;
            % Step 5: Objectives (mixed)    
            L      = 2; % number of segments
            G      = 1-[ones(N,1) cumprod(sin(pi/2*THETA),2)] .* [cos(pi/2*THETA) ones(N,1)];
            G(:,1) = THETA(:,1) - cos(2*pi*L*THETA(:,1)+pi/2)/2/L/pi;
            PopObj = G .* repmat((1+T),1,M) ;
            % Step 6: Constraints
            PopCon(:,1) = Sx(:,1) + h - 1; 
            PopCon(:,2) = -(Sx(:,1) + h - 3/4);
            Population  = SOLUTION(varargin{1},PopObj,PopCon,varargin{2:end});
            obj.FE      = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,obj.M);
            R = 1 - R./repmat(sqrt(sum(R.^2,2)),1,obj.M);
            V = asin(sqrt(sum((1-R(:,2:end)).^2,2)));
            R(:,1) = (2/pi).*(V+1/8*sin(8*V));
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                R = obj.GetOptimum(100);
            elseif obj.M == 3
                a = linspace(0,pi/2,10)';
                R = {1-sin(a)*cos(a'),1-sin(a)*sin(a'),1-cos(a)*ones(size(a'))};
                V = asin(sqrt((1-R{2}).^2+(1-R{3}).^2));
                R{1} = (2/pi).*(V+1/8*sin(8*V));
            else
                R = [];
            end
        end
    end
end