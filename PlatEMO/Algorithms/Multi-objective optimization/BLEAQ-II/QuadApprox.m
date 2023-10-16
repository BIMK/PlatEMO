function approxmodel = QuadApprox(y,X)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    warning off all;
    dimensions  = size(X,2);
    datasetSize = size(X,1);
    if (datasetSize<dimensions)
        error('Cannot compute model as the datasetSize is smaller than dimensions.');
    end

    modelConsidered = {'linear','purequadratic','quadratic'};
    numModel = length(modelConsidered);
    XX       = cell(1,numModel);
    output   = cell(1,numModel);
    mse      = zeros(1,numModel);
    bic      = zeros(1,numModel);

    for i = 1 : numModel
        model = modelConsidered{i};
        XX{i} = x2fx(X, model);
        output{i}.beta = XX{i}\y; 
        mse(i) = real(sum((y-XX{i}*output{i}.beta).^2)/length(y));
        bic(i) = size(XX{i},2)*log(size(XX{i},1))/size(XX{i},1) + 2*log(mse(i)); 
    end

    [~,index] = min(bic);

    modelSelected = modelConsidered(index);
    constant      = output{index}.beta(1);
    linear        = output{index}.beta(2:2+dimensions-1);
    sqmatrix      = zeros(dimensions,dimensions);

    if strcmp(modelSelected,'purequadratic')
        diagonal = output{i}.beta(end-dimensions+1:end);
        for i = 1 : dimensions
            sqmatrix(i,i) = diagonal(i);
        end
    end

    if strcmp(modelSelected,'quadratic')
        cross    = output{i}.beta(2+dimensions:end-dimensions);
        diagonal = output{i}.beta(end-dimensions+1:end);
        k        = 0;
        for i = 1 : dimensions
            sqmatrix(i,i) = diagonal(i);
            for j = i+1 : dimensions
                k = k + 1;
                sqmatrix(i,j) = cross(k)/2;
                sqmatrix(j,i) = cross(k)/2;
            end
        end
    end

    stdXX           = std(XX{index});
    stdy            = std(y);
    stdXX(stdXX==0) = -realmin;
    stdy(stdy==0)   = -realmin;
    XXNorm          = (XX{index}-ones(size(XX{index},1),1)*mean(XX{index}))./(ones(size(XX{index},1),1)*stdXX);
    yNorm           = (y-mean(y))./stdy;
    betaNorm        = XXNorm\yNorm;
    mseNorm         = sum((yNorm-XXNorm*betaNorm).^2)/length(yNorm);

    approxmodel = struct('model',modelSelected,'constant',real(constant),...
                    'linear',real(linear),'sqmatrix',real(sqmatrix),...
                    'mse',mse(index),'mseNorm',mseNorm,'bic',bic(index));

    if isnan(approxmodel.mse)
        error('The code has ended up with NAN values. One of the reasons could be duplicate rows in the input matrix that makes the system under-defined and a solution cannot be found using least-squares.')
    end
end