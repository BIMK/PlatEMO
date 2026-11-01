classdef MASSRFeatures
%MASSRFeatures - Population observations and auxiliary guidance.
%
%   All methods are static. Observation columns contain bound-normalized
%   decisions and population-normalized objectives. Guidance targets
%   supervise representation learning and are not part of the deployed
%   state.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods(Static)
        function X = observation(Population,Problem)
        %observation - Build the normalized population observation.
        %
        %   X = MASSRFeatures.observation(Population,Problem) returns a
        %   (D+M)-by-N matrix. Decisions use problem bounds and objectives
        %   use the current population minimum and maximum.

            Dec           = double(Population.decs);
            lower         = reshape(Problem.lower,1,[]);
            upper         = reshape(Problem.upper,1,[]);
            span          = upper - lower;
            span(span<=0) = 1;
            Dec           = min(max((Dec-lower)./span,0),1);
            X             = [Dec,MASSRFeatures.normalize(double(Population.objs))]';
        end

        function G = guidance(Population,LastPopulation,Problem,improvement,stagnation,replacementRate)
        %guidance - Extract the 12 auxiliary supervision targets.
        %
        %   G = MASSRFeatures.guidance(...) returns 12 progress,
        %   improvement, quality, diversity, spread, movement, and
        %   replacement indicators in the order expected by the head masks.
        %
        %   improvement is a signed relative change; stagnation and
        %   replacementRate are supplied by the caller. Targets are clipped
        %   to [0,1]. Initial movement indicators are zero when
        %   LastPopulation is empty.

            X   = MASSRFeatures.observation(Population,Problem)';
            Dec = X(:,1:Problem.D);
            Obj = X(:,Problem.D+1:end);
            raw = double(Population.objs);

            % Progress, improvement, and high-quality solution ratio
            progress    = Problem.FE/Problem.maxFE;
            improvement = 0.5 + 0.5*tanh(5*improvement);
            ratio       = mean(NDSort(raw,1)==1);

            % Diversity and spread in objective and decision spaces
            objDiversity = pairDistance(Obj);
            uniformity   = distribution(Obj);
            objSpread    = mean(std(Obj,0,1));
            decDiversity = pairDistance(Dec);
            decSpread    = mean(std(Dec,0,1));

            % Movement between consecutive populations and objective bounds
            movement      = 0;
            frontMovement = 0;
            if ~isempty(LastPopulation)
                last          = MASSRFeatures.observation(LastPopulation,Problem)';
                L             = min(size(Dec,1),size(last,1));
                movement      = mean(sqrt(sum((Dec(1:L,:)-last(1:L,1:Problem.D)).^2,2)))/sqrt(Problem.D);
                oldObj        = double(LastPopulation.objs);
                lo            = min(raw,[],1);
                oldLo         = min(oldObj,[],1);
                hi            = max(raw,[],1);
                oldHi         = max(oldObj,[],1);
                frontMovement = 0.5*mean(abs(lo-oldLo)./(abs(lo)+abs(oldLo)+1e-12)) + ...
                    0.5*mean(abs(hi-oldHi)./(abs(hi)+abs(oldHi)+1e-12));
            end

            % Preserve the guidance order used by all head masks
            G               = [progress;improvement;ratio;objDiversity;uniformity;objSpread; ...
                decDiversity;decSpread;movement;stagnation;replacementRate;frontMovement];
            G(~isfinite(G)) = 0;
            G               = min(max(G,0),1);
        end

        function mask = mask(testbed)
        %mask - Return the role-specific guidance masks.
        %
        %   M = MASSRFeatures.mask(testbed) returns a 12-by-K binary matrix.
        %   Each column selects the guidance components used to supervise
        %   one readout head.

            switch upper(testbed)
                case 'MADE'
                    indices = {[1,2,3,7,8,9],[1,2,4,6,9,11],[1,2,7,10,11,12]};
                case 'MAMOEAD'
                    indices = {[1,2,3,4,5,6],[1,4,5,7,8,9],[1,2,7,8,9,11],[1,2,7,10,11,12]};
                case 'CMODRL'
                    indices = {[1,2,3,4,5,10],[1,2,7,8,9,11]};
                case 'SPARSE'
                    indices = {[1,2,3,7,8,10],[1,3,7,8,9,10],[1,2,7,8,9,11]};
            end
            mask = zeros(12,numel(indices));
            for k = 1 : numel(indices)
                mask(indices{k},k) = 1;
            end
        end

        function value = change(before,after)
        %change - Compute relative improvement in a minimized value.
        %
        %   The result is (before-after)/max(abs(before),1e-12). Positive
        %   values indicate improvement.

            value = (before-after)/max(abs(before),1e-12);
        end

        function D = distance(A,B)
        %distance - Compute Euclidean distances between rows.
        %
        %   D = MASSRFeatures.distance(A,B) has size size(A,1)-by-size(B,1).
        %   Nonnegative clipping removes roundoff in squared distances.

            D = sqrt(max(sum(A.^2,2)+sum(B.^2,2)'-2*(A*B'),0));
        end

        function Y = normalize(X)
        %normalize - Normalize each feature to the unit interval.
        %
        %   Constant features map to zero. Other columns use their observed
        %   minimum and maximum.

            lower            = min(X,[],1);
            span             = max(X,[],1) - lower;
            span(span<1e-12) = 1;
            Y                = min(max((X-lower)./span,0),1);
        end
    end
end

function value = pairDistance(X)
%pairDistance - Measure mean pairwise population diversity.
%
%   The mean distance over distinct pairs is divided by the square root of
%   the feature dimension.

    N = size(X,1);
    if N < 2
        value = 0;
        return;
    end
    D     = MASSRFeatures.distance(X,X);
    value = sum(D(triu(true(N),1)))/(N*(N-1)/2)/sqrt(size(X,2));
end

function value = distribution(X)
%distribution - Measure nearest-neighbor distance uniformity.
%
%   Use the coefficient of variation of nearest-neighbor distances. Larger
%   values indicate more uniform spacing.

    N = size(X,1);
    if N < 3
        value = 0;
        return;
    end
    D            = MASSRFeatures.distance(X,X);
    D(1:N+1:end) = inf;
    nearest      = min(D,[],2);
    value        = 1/(1+std(nearest)/(mean(nearest)+1e-12));
end