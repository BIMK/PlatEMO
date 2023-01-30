classdef Sparse_SR < PROBLEM
% <multi> <real> <large/none> <expensive/none> <sparse/none>
% The sparse signal reconstruction problem
% lenSig   --- 1024 --- Length of signal 
% lenObs   ---  480 --- Number of observations
% sparsity ---  260 --- Sparsity  
% sigma    ---    0 --- Noise level

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

    properties(Access = private)
        A;          % Measurement matrix
        b;          % Observations
        x_true;		% True signal
        Total;		
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            [lenSig,lenObs,sparsity,sigma] = obj.ParameterSet(1024,480,260,0);
            N = lenSig;
            M = lenObs;
            K = sparsity;
            fileName = fullfile(fileparts(mfilename('fullpath')),sprintf('Dataset_SR-N%d-M%d-K%d-sigma%.2f.mat',N,M,K,sigma));
            if exist(fileName,'file') == 2
                load(fileName);
            else
                % Generate dataset
                [A, b, x_true] = inst_gen(N, M, K, sigma);
               save(fileName,'A','b','x_true','N','M','K','sigma');
            end
            [obj.A,obj.b,obj.x_true] = deal(A,b,x_true);
            % Parameter setting
            obj.M = 2;
            obj.D = N;
            x = obj.x_true(obj.x_true~=0);           
            obj.lower    = floor(repmat(mean(x)-3*std(x),1,obj.D));
            obj.upper    = ceil(ones(1,obj.D)*(mean(x)+3*std(x)));
            obj.encoding = ones(1,obj.D);

            x = zeros(1, N);
            obj.Total = norm(obj.A*x'-obj.b);
        end
        %% Generate initial solutions
        function Population = Initialization(obj,N)
            if nargin < 2
                N = obj.N;
            end
            e = randi(obj.D,1,N);
            PopDec = zeros(N,obj.D);
            for i = 1 : N
                integer        = e(i);          % sparse support
                [~,spar_b]     = sort(rand(1,obj.D));
                sparse_support = spar_b(1:integer);
                SS = randn(integer,1); 
                SS = SS./norm(SS);
                PopDec(i,sparse_support) = SS;	% initialize sparse signal
            end
            Population = obj.Evaluation(PopDec);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopDec = logical(PopDec);
            PopObj = zeros(size(PopDec,1),2);
            [N,D]  = size(PopDec);
            PopObj(:,1) = sum(PopDec~=0,2)/D;
            for i = 1:N
                PopObj(i,2) = norm(obj.A*PopDec(i,:)'-obj.b) / obj.Total;
            end
        end
        %% Display a population in the objective space
        function DrawObj(obj,Population)
            Draw(Population.objs,{'Sparsity of reconstructed signal','Loss of reconstructed signal',[]});
        end
    end
end

function [A,b,x_true] = inst_gen(N,M,K,sigma)
    x_true = zeros(N,1);
    q = randperm(N);
    x_true(q(1:K)) = 2*randn(K,1);

    err = sigma * randn(M,1);

    A = randn(M,N);
    A = orth(A')';
    b = A*x_true;

    normA = norm(A);
    A = A/normA;
    b = b/normA + err;
end