classdef MASSREncoder < handle
%MASSREncoder - Self-supervised population encoder.
%
%   One encoder is trained per optimization run from early population
%   observations. Inputs are (D+M)-by-N matrices, with individuals in
%   columns. Each independent readout head returns a 16-dimensional state.
%
%   collect stores observations and guidance targets. fit trains the
%   backbone and readouts once, then discards the training heads and
%   samples. encode uses only population observations after training.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = private)
        Model;            % Backbone and state/training head parameters
        Ready = false;    % Whether one-time representation training is complete
    end
    properties(Access = private)
        Mask;                 % 12-by-K role-specific guidance masks
        Observations = {};    % Collected C-by-N population matrices
        Guidance     = [];    % 12-by-B targets for the collected observations
    end

    methods
        function obj = MASSREncoder(observationDimension,testbed)
        %MASSREncoder - Create the backbone and independent readout heads.
        %
        %   obj = MASSREncoder(C,testbed) creates a C-input encoder. testbed
        %   is MADE, MAMOEAD, CMODRL, or SPARSE and determines the number of
        %   heads and their guidance masks.

            obj.Mask  = MASSRFeatures.mask(testbed);
            obj.Model = MASSREncoder.initialize(observationDimension,size(obj.Mask,2));
        end

        function collect(obj,X,G)
        %collect - Store one population observation and its guidance.
        %
        %   obj.collect(X,G) appends a C-by-N observation and a 12-by-1
        %   guidance target. Calls belong to the early random-configuration
        %   phase; collection does not evaluate new solutions.

            obj.Observations{end+1} = X;
            obj.Guidance(:,end+1)   = G;
        end

        function fit(obj)
        %fit - Train once and retain the inference module.
        %
        %   obj.fit() optimizes the backbone and all heads for 40 epochs
        %   using two augmented views and the original-view guidance
        %   targets. Projection heads, guidance predictors, and stored
        %   samples are then released.

            X         = obj.Observations;
            G         = obj.Guidance;
            optimizer = [];

            % Optimize the decaying guidance and contrastive objectives
            for epoch = 1 : 40
                alpha = exp(-3*(epoch-1)/39);
                order = randperm(numel(X));
                first = 1;
                while first <= numel(order)
                    count = min(32,numel(order)-first+1);
                    if numel(order) - first + 1 - count == 1
                        count = count - 1;
                    end
                    ids                   = order(first:first+count-1);
                    views.raw             = X(ids);
                    views.first           = cellfun(@augment,X(ids),'UniformOutput',false);
                    views.second          = cellfun(@augment,X(ids),'UniformOutput',false);
                    [~,gradient]          = MASSREncoder.objective(obj.Model,views,G(:,ids),alpha,obj.Mask);
                    [obj.Model,optimizer] = MASSRNN.adam(obj.Model,gradient,optimizer,1e-3);
                    first                 = first + count;
                end
            end

            % Retain only the frozen backbone and readout heads
            obj.Model        = rmfield(obj.Model,{'projection','guidance'});
            obj.Ready        = true;
            obj.Observations = {};
            obj.Guidance     = [];
        end

        function state = encode(obj,X)
        %encode - Return the agent-specific population states.
        %
        %   S = obj.encode(X) returns a 16-by-K matrix from one C-by-N
        %   observation. Column k is the state for agent k. Guidance targets
        %   and rewards are not inputs to this method.

            tokens = MASSREncoder.tokens(obj.Model,X);
            state  = zeros(16,numel(obj.Model.head));
            for k = 1 : numel(obj.Model.head)
                state(:,k) = readout(obj.Model.head{k},tokens);
            end
        end
    end

    methods(Static)
        function P = initialize(dimension,heads)
        %initialize - Initialize the shared backbone and K sets of heads.
        %
        %   P = MASSREncoder.initialize(C,K) creates 32-dimensional
        %   embeddings and tokens, K attention readouts, 16-dimensional
        %   projection outputs, and 12-dimensional guidance predictors.

            P.phi        = MASSRNN.initialize([dimension,32]);
            P.score      = MASSRNN.initialize([32,1]);
            P.psi        = MASSRNN.initialize([32,32,32]);
            P.head       = cell(1,heads);
            P.projection = cell(1,heads);
            P.guidance   = cell(1,heads);
            for k = 1 : heads
                P.head{k}       = struct('q',(2*rand(32,1)-1)*sqrt(6/33), ...
                    'W',(2*rand(16,32)-1)*sqrt(6/48),'b',zeros(16,1));
                P.projection{k} = MASSRNN.initialize([16,16,16]);
                P.guidance{k}   = MASSRNN.initialize([16,12]);
            end
        end

        function [Z,cache] = tokens(P,X)
        %tokens - Construct 12 representative population tokens.
        %
        %   [Z,cache] = MASSREncoder.tokens(P,X) returns a 32-by-12 token
        %   matrix. Attention weights rescale individual embeddings;
        %   distance scores in normalized observation space select the
        %   representatives.

            % Embed individuals and compute population attention
            [H,cache.phi]   = MASSRNN.forward(P.phi,double(X),true);
            [e,cache.score] = MASSRNN.forward(P.score,H);
            omega           = softmaxRow(e);

            % Select representatives in normalized observation space
            distances                            = MASSRFeatures.distance(X',X');
            distances(1:size(distances,1)+1:end) = 0;
            key                                  = [-sum(distances,2),double(X')];
            [~,order]                            = sortrows(key);
            index                                = order(1:12);

            % Project the selected attention-weighted embeddings
            weighted      = H(:,index).*omega(index);
            [Z,cache.psi] = MASSRNN.forward(P.psi,weighted);
            cache.H       = H;
            cache.omega   = omega;
            cache.index   = index;
        end

        function [loss,G] = objective(P,views,target,alpha,mask)
        %objective - Evaluate both self-supervised losses and gradients.
        %
        %   [loss,G] = MASSREncoder.objective(P,views,target,alpha,mask)
        %   sums the head-specific contrastive and masked guidance losses.
        %   views contains raw, first, and second cell arrays; target has 12
        %   rows.

            B     = numel(views.raw);
            input = [views.raw,views.first,views.second];
            T     = cell(1,3*B);
            cache = cell(size(T));
            dT    = cell(size(T));

            % Reuse backbone tokens across all head-specific objectives
            for i = 1 : 3*B
                [T{i},cache{i}] = MASSREncoder.tokens(P,input{i});
                dT{i}           = zeros(size(T{i}));
            end
            G    = MASSRNN.zeros(P);
            loss = 0;
            for k = 1 : numel(P.head)
                S         = zeros(16,3*B);
                headCache = cell(1,3*B);
                for i = 1 : 3*B
                    [S(:,i),headCache{i}] = readout(P.head{k},T{i});
                end

                % Contrastive supervision uses the two augmented views
                [projection,projectCache]        = MASSRNN.forward(P.projection{k},S(:,B+1:end));
                [contrast,dProjection]           = MASSREncoder.contrastive(projection,0.2);
                [G.projection{k},dContrastState] = MASSRNN.backward( ...
                    P.projection{k},projectCache,(1-alpha)*dProjection);

                % Masked guidance supervision uses the original view
                [prediction,guideCache]     = MASSRNN.forward(P.guidance{k},S(:,1:B));
                residual                    = (prediction-target).*mask(:,k);
                divisor                     = B*max(sum(mask(:,k)),1);
                guide                       = sum(residual(:).^2)/divisor;
                [G.guidance{k},dGuideState] = MASSRNN.backward( ...
                    P.guidance{k},guideCache,2*alpha*residual/divisor);
                dState                      = [dGuideState,dContrastState];
                for i = 1 : 3*B
                    [gHead,gTokens] = readoutBackward(P.head{k},headCache{i},dState(:,i));
                    G.head{k}       = MASSRNN.add(G.head{k},gHead);
                    dT{i}           = dT{i} + gTokens;
                end
                loss = loss + (1-alpha)*contrast + alpha*guide;
            end
            for i = 1 : 3*B
                % Backpropagate the accumulated token gradients through the backbone
                gBackbone = tokensBackward(P,cache{i},dT{i});
                G.phi     = MASSRNN.add(G.phi,gBackbone.phi);
                G.score   = MASSRNN.add(G.score,gBackbone.score);
                G.psi     = MASSRNN.add(G.psi,gBackbone.psi);
            end
        end

        function [loss,gradient] = contrastive(P,tau)
        %contrastive - Evaluate paired-view InfoNCE with cosine similarity.
        %
        %   [loss,G] = MASSREncoder.contrastive(P,tau) treats columns of P
        %   as 2B projected views. The first B columns are paired with the
        %   final B columns. Self-comparisons are excluded and each anchor
        %   has 2B-2 negative views.

            count                 = size(P,2);
            norms                 = sqrt(sum(P.^2,1));
            scale                 = max(norms,1e-12);
            U                     = P./scale;
            logits                = (U'*U)/tau;
            logits(1:count+1:end) = -inf;
            positive              = [count/2+1:count,1:count/2];
            rows                  = (1:count)';
            ids                   = sub2ind([count,count],rows,positive');
            maximum               = max(logits,[],2);
            exponential           = exp(logits-maximum);
            denominator           = sum(exponential,2);
            loss                  = mean(maximum+log(denominator)-logits(ids));
            dLogits               = exponential./denominator;
            dLogits(ids)          = dLogits(ids) - 1;
            dLogits               = dLogits/count;
            dU                    = U*(dLogits+dLogits')/tau;
            gradient              = (dU-U.*sum(U.*dU,1))./scale;
            small                 = norms < 1e-12;
            gradient(:,small)     = dU(:,small)/1e-12;
        end
    end
end

function X = augment(X)
%augment - Subsample, perturb, and mask a population view.
%
%   Keep 85% of individuals, with at least 12 representatives. Add
%   feature-scaled Gaussian noise, mask 5% of entries, and clip to [0,1].

    N     = size(X,2);
    X     = double(X(:,randperm(N,max(12,ceil(0.85*N)))));
    sigma = 0.03*(std(X,0,2)+1e-6);
    X     = X + sigma.*randn(size(X));
    X     = X.*(rand(size(X))>=0.05);
    X     = min(max(X,0),1);
end

function [state,cache] = readout(head,T)
%readout - Pool shared tokens for one agent.
%
%   The learned query defines attention over tokens. An affine output map
%   converts the pooled 32-dimensional vector to a 16-dimensional state.

    beta   = softmaxRow(head.q'*T);
    pooled = T*beta';
    state  = head.W*pooled + head.b;
    cache  = struct('T',T,'beta',beta,'pooled',pooled);
end

function [G,dT] = readoutBackward(head,cache,ds)
%readoutBackward - Differentiate one attention readout.
%
%   Return gradients for the query and affine map, together with the
%   gradient propagated to the shared tokens.

    G.W   = ds*cache.pooled';
    G.b   = ds;
    dp    = head.W'*ds;
    dbeta = dp'*cache.T;
    de    = cache.beta.*(dbeta-sum(dbeta.*cache.beta));
    G.q   = cache.T*de';
    dT    = dp*cache.beta + head.q*de;
end

function G = tokensBackward(P,cache,dT)
%tokensBackward - Backpropagate through token construction.
%
%   Gradients pass through the selected embeddings, population softmax
%   scores, and individual embedding network. Selection indices remain fixed
%   for the current forward pass.

    [G.psi,dWeighted]   = MASSRNN.backward(P.psi,cache.psi,dT);
    dH                  = zeros(size(cache.H));
    dH(:,cache.index)   = dWeighted.*cache.omega(cache.index);
    dOmega              = zeros(size(cache.omega));
    dOmega(cache.index) = sum(dWeighted.*cache.H(:,cache.index),1);
    de                  = cache.omega.*(dOmega-sum(dOmega.*cache.omega));
    [G.score,dScore]    = MASSRNN.backward(P.score,cache.score,de);
    [G.phi,~]           = MASSRNN.backward(P.phi,cache.phi,dH+dScore);
end

function p = softmaxRow(x)
%softmaxRow - Compute a stable softmax across each row.

    p = exp(x-max(x,[],2));
    p = p./sum(p,2);
end