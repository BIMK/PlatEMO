% SAVECLASS - Saves the SVM in memory to the specified file.
%
% Syntax: saveclass(fname)
%
% Version 3.22e -- Comments to diehl@alumni.cmu.edu
%

function saveclass(fname)

% define global variables 
global a;                     % alpha coefficients
global b;                     % bias
global C;                     % regularization parameters
global deps;                  % jitter factor in kernel matrix
global g;                     % partial derivatives of cost function w.r.t. alpha coefficients
global ind;                   % cell array containing indices of margin, error, reserve and unlearned vectors
global max_reserve_vectors;   % maximum number of reserve vectors stored
global Q;                     % extended kernel matrix for all vectors
global Rs;                    % inverse of extended kernel matrix for margin vectors   
global scale;                 % kernel scale
global type;                  % kernel type
global uind;                  % user-defined example indices
global X;                     % matrix of margin, error, reserve and unlearned vectors stored columnwise
global y;                     % column vector of class labels (-1/+1) for margin, error, reserve and unlearned vectors

% save SVM state
cmd = sprintf('save %s a b C deps g ind max_reserve_vectors Q Rs scale type uind X y;',fname);
eval(cmd);
