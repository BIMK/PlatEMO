function [validator, estimator] = generateFunctions(M,alpha)
% Generate functions

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Tomoaki Takagi

    %% Generate validator
    if M > 4 || isinf(alpha)
        validator = @(W,X) W;
    elseif M == 2
        validator = @clusterValidator;
    else
        validator = @alphaShapeValidator;
    end

    %% Generate estimator
    if ~isempty(ver('nnet'))
        estimator = @(W,X,Y) W.*sim(newrbe(X',Y'),W')';
    else
        I = ones(1,M);
        estimator = @(W,X,Y) W.*predictor(W, ...
            dacefit(X,Y,'regpoly0','corrgauss',I,1e-3*I,1e3*I));
    end

    %% Validator
    i = 1:M-1;
    A = [0*i;
        tril(repmat(sqrt(1./i./(i+1)),M-1,1),-1)+diag(sqrt((i+1)./i))];

    function W = clusterValidator(W,X)
        XA = sort(X*A);
        T  = clusterdata(XA,'Criterion','distance','Cutoff',alpha);
        XA = XA(setxor(uni(T),uni(T,'last')));
        W  = W(mod(discretize(W*A,XA),2)==1,:);
    end

    function W = alphaShapeValidator(W,X)
        shp = alphaShape(X*A,alpha);
        if alpha > min(alphaSpectrum(shp))
            W = W(inShape(shp,W*A),:);
        end
    end
end

function ia = uni(varargin)
    [~,ia] = unique(varargin{:});
end 