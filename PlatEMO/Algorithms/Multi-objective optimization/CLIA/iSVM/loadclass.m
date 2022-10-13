% LOADCLASS - Loads the SVM from the specified file into memory.
%
% Syntax: loadclass(fname)
%
% Version 3.22e -- Comments to diehl@alumni.cmu.edu
%

function loadclass(fname)

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

% load SVM state
flag = 0;
cmd = sprintf('load %s a b C deps g ind max_reserve_vectors Q Rs scale type uind X y;',fname);
eval(cmd,'flag = 1;');
if (flag)
   disp('LOADCLASS: Unable to load the SVM!');
end;
