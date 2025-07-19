function OffspringLearn = DA_TwoLayer(SourceTemp,SourceGmax,Target,Problem,D)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    % input should be saved in column
    Noisy = SourceTemp'; % noisy signals
    
    Clean = SourceGmax'; % clean signals
    % clean signals
    % add constant feature to the input, x = [x;1]
    [d,n]    = size(Noisy);
    NoisyAdd = [Noisy;ones(1,n)];
    CleanAdd = [Clean;ones(1,n)];
    
    % calculate the mapping
    P = CleanAdd*NoisyAdd';
    Q = NoisyAdd*NoisyAdd';
    lambda       = 1e-5;
    reg          = lambda*eye(d+1);
    reg(end,end) = 0;
    W1           = P/(Q+reg);
    
    % the second layer
    NoisyAdd2 = tanh(W1*NoisyAdd);
    P  = CleanAdd*NoisyAdd2';
    Q  = NoisyAdd2*NoisyAdd2';
    W2 = P/(Q+reg);
    % delete the bias in the mapping, i.e., the b in M=[M,b]
    b       = size(W1,1);
    W1(b,:) = [];
    W1(:,b) = [];
    W2(b,:) = [];
    W2(:,b) = [];
    % apply the mapping to the test signals
    OffspringDec = (W2*tanh(W1*Target'))';
    % boundary check
    lower = [0,-ones(1,D-1)];
    upper = [1, ones(1,D-1)];
    Lower = repmat(lower,size(OffspringDec,1),1);
    Upper = repmat(upper,size(OffspringDec,1),1);
    OffspringDec   = max(min(OffspringDec,Upper),Lower); % N*D
    OffspringLearn = Problem.Evaluation(OffspringDec);
end