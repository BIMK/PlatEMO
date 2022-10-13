% UPDATERQ - Updates Rs and Q accordingly when adding or removing a 
%            margin vector.  Note: Rs and Q are defined as global
%            variables.    
%
% Syntax: updateRQ(beta,gamma,indc) 
%         (for adding a margin vector)   
%         
%         updateRQ(indc)
%         (for removing a margin vector)
%
%   beta: parameter sensitivities associated with the example indc 
%  gamma: margin sensitivity associated with the example indc
%   indc: example/matrix index
%
% Version 3.22e -- Comments to diehl@alumni.cmu.edu
%

function updateRQ(varargin)

% flags for example state
MARGIN    = 1;
ERROR     = 2;
RESERVE   = 3;
UNLABELED = 4;

% define global variables
global C; 		% regularization parameters
global ind;    % cell array containing indices of margin, error, reserve and unlearned vectors
global deps;   % jitter factor in kernel matrix
global Q;      % extended kernel matrix for all vectors
global Rs;     % inverse of extended kernel matrix for margin vectors   
global scale;  % kernel scale
global type;   % kernel type
global X;      % matrix of margin, error and reserve vectors stored columnwise
global y;      % column vector of class labels (-1/+1) for margin, error and reserve vectors

if (nargin == 3)
   beta  = varargin{1};
   gamma = varargin{2};
   indc = varargin{3};
   expand = 1;   
elseif (nargin == 1)
   indc = varargin{1};
   expand = 0;   
else
   error('updateRQ: Incorrect number of parameters');
end;

rows = size(Rs,1);
if (expand)
   
   if (gamma < deps)
      gamma = deps;
   end;
   
   if (rows > 1)
      Rs = [Rs,zeros(rows,1);zeros(1,rows+1)] + [beta;1]*[beta',1]/gamma;
   else
      Rs = [-(kernel(X(:,indc),X(:,indc),type,scale)+deps) y(indc) ; y(indc) 0];
   end;
   Q = [Q ; (y(indc)*y').*kernel(X(:,indc),X,type,scale)];
   Q(rows+1,indc) = Q(rows+1,indc) + deps;
   
else
   
   if (rows > 2)
      stripped = [1:indc-1 indc+1:size(Rs,1)]; 
      Rs = Rs(stripped,stripped)-Rs(stripped,indc)*Rs(indc,stripped)/Rs(indc,indc);
   else
      Rs = Inf;
   end;
   Q(indc,:) = [];
     
end;

