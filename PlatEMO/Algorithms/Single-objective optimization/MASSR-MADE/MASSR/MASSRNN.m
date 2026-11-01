classdef MASSRNN
%MASSRNN - Dense networks, backpropagation, and Adam.
%
%   All methods are static and operate on numeric arrays, structures, or
%   cell arrays. Networks store weights and biases in net.W and net.b.
%   Samples are arranged in columns; hidden layers use ReLU.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods(Static)
        function net = initialize(widths)
        %initialize - Initialize a dense network with zero biases.
        %
        %   net = MASSRNN.initialize(widths) creates one affine layer per
        %   adjacent pair of widths, using Glorot uniform weights.

            net.W = cell(1,numel(widths)-1);
            net.b = cell(size(net.W));
            for i = 1 : numel(net.W)
                limit    = sqrt(6/(widths(i)+widths(i+1)));
                net.W{i} = (2*rand(widths(i+1),widths(i))-1)*limit;
                net.b{i} = zeros(widths(i+1),1);
            end
        end

        function [Y,cache] = forward(net,X,lastRelu)
        %forward - Evaluate a network and cache layer inputs.
        %
        %   [Y,cache] = MASSRNN.forward(net,X,lastRelu) uses samples in
        %   columns. Hidden layers use ReLU. lastRelu defaults to false and
        %   controls the final-layer activation.

            if nargin < 3
                lastRelu = false;
            end
            cache.X    = cell(size(net.W));
            cache.Z    = cell(size(net.W));
            cache.relu = [true(1,numel(net.W)-1),lastRelu];
            Y          = X;
            for i = 1 : numel(net.W)
                cache.X{i} = Y;
                Y          = net.W{i}*Y + net.b{i};
                cache.Z{i} = Y;
                if cache.relu(i)
                    Y = max(Y,0);
                end
            end
        end

        function [gradient,dX] = backward(net,cache,dY)
        %backward - Differentiate a cached network evaluation.
        %
        %   [G,dX] = MASSRNN.backward(net,cache,dY) returns the parameter
        %   gradients and the input gradient for the supplied output
        %   derivative.

            gradient = MASSRNN.zeros(net);
            for i = numel(net.W) : -1 : 1
                if cache.relu(i)
                    dY = dY.*(cache.Z{i}>0);
                end
                gradient.W{i} = dY*cache.X{i}';
                gradient.b{i} = sum(dY,2);
                dY            = net.W{i}'*dY;
            end
            dX = dY;
        end

        function Y = zeros(X)
        %zeros - Create a zero-valued parameter tree.
        %
        %   Structures and cell arrays retain their organization; numeric
        %   leaves are replaced by zero arrays of the same size.

            if isstruct(X)
                Y     = struct();
                names = fieldnames(X);
                for i = 1 : numel(names)
                    Y.(names{i}) = MASSRNN.zeros(X.(names{i}));
                end
            elseif iscell(X)
                Y = cell(size(X));
                for i = 1 : numel(X)
                    Y{i} = MASSRNN.zeros(X{i});
                end
            else
                Y = zeros(size(X));
            end
        end

        function C = add(A,B)
        %add - Add two identically structured parameter trees.

            if isstruct(A)
                C     = A;
                names = fieldnames(A);
                for i = 1 : numel(names)
                    C.(names{i}) = MASSRNN.add(A.(names{i}),B.(names{i}));
                end
            elseif iscell(A)
                C = cell(size(A));
                for i = 1 : numel(A)
                    C{i} = MASSRNN.add(A{i},B{i});
                end
            else
                C = A + B;
            end
        end

        function Y = scale(X,factor)
        %scale - Multiply every numeric leaf by a scalar.

            if isstruct(X)
                Y     = X;
                names = fieldnames(X);
                for i = 1 : numel(names)
                    Y.(names{i}) = MASSRNN.scale(X.(names{i}),factor);
                end
            elseif iscell(X)
                Y = cell(size(X));
                for i = 1 : numel(X)
                    Y{i} = MASSRNN.scale(X{i},factor);
                end
            else
                Y = factor*X;
            end
        end

        function value = squaredNorm(X)
        %squaredNorm - Sum squared values across all numeric leaves.

            if isstruct(X)
                value = 0;
                names = fieldnames(X);
                for i = 1 : numel(names)
                    value = value + MASSRNN.squaredNorm(X.(names{i}));
                end
            elseif iscell(X)
                value = 0;
                for i = 1 : numel(X)
                    value = value + MASSRNN.squaredNorm(X{i});
                end
            else
                value = sum(X(:).^2);
            end
        end

        function [P,state] = adam(P,G,state,rate,maxNorm)
        %adam - Apply Adam to a nested parameter tree.
        %
        %   [P,state] = MASSRNN.adam(P,G,state,rate,maxNorm) updates
        %   parameters and optimizer moments. An empty state initializes the
        %   moments. maxNorm defaults to infinity and optionally clips the
        %   global gradient norm.

            if nargin < 5
                maxNorm = inf;
            end
            normG = sqrt(MASSRNN.squaredNorm(G));
            if normG > maxNorm
                G = MASSRNN.scale(G,maxNorm/normG);
            end
            if isempty(state)
                state   = struct();
                state.m = MASSRNN.zeros(P);
                state.v = MASSRNN.zeros(P);
                state.t = 0;
            end
            state.t             = state.t + 1;
            [P,state.m,state.v] = updateTree(P,G,state.m,state.v,state.t,rate);
        end
    end
end

function [P,m,v] = updateTree(P,G,m,v,t,rate)
%updateTree - Update Adam moments and parameters recursively.
%
%   Numeric leaves use first- and second-moment decay rates of 0.9 and 0.999
%   with bias correction at update t.

    if isstruct(P)
        names = fieldnames(P);
        for i = 1 : numel(names)
            key                       = names{i};
            [P.(key),m.(key),v.(key)] = updateTree(P.(key),G.(key),m.(key),v.(key),t,rate);
        end
    elseif iscell(P)
        for i = 1 : numel(P)
            [P{i},m{i},v{i}] = updateTree(P{i},G{i},m{i},v{i},t,rate);
        end
    else
        m = 0.9*m + 0.1*G;
        v = 0.999*v + 0.001*G.^2;
        P = P - rate*(m/(1-0.9^t))./(sqrt(v/(1-0.999^t))+1e-8);
    end
end