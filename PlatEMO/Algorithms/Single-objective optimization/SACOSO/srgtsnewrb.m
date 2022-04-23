function net = srgtsnewrb(p,t,goal,spread,mn,df,displayLevel)
%Function srgtsnewrb designs a radial basis network. This is just a modificated
%version of the "newrb" MATLAB function. The only modification consist in add
%display levels to the training process to match with the SurrogatesToolBox.
%
%  Synopsis
%
%    net = srgtsnewrb
%    [net,tr] = srgtsnewrb(P,T,GOAL,SPREAD,MN,DF)
%
%  Description
%
%    Radial basis networks can be used to approximate
%    functions.  NEWRB adds neurons to the hidden
%    layer of a radial basis network until it meets
%    the specified mean squared error goal.
%
%   NET = NEWRB creates a new network with a dialog box.
%
%   NEWRB(PR,T,GOAL,SPREAD,MN,DF) takes these arguments,
%     P      - RxQ matrix of Q input vectors.
%     T      - SxQ matrix of Q target class vectors.
%     GOAL   - Mean squared error goal, default = 0.0.
%     SPREAD - Spread of radial basis functions, default = 1.0.
%     MN     - Maximum number of neurons, default is Q.
%     DF     - Number of neurons to add between displays, default = 25.
%   and returns a new radial basis network.
%
%   The larger that SPREAD is the smoother the function approximation
%   will be.  Too large a spread means a lot of neurons will be
%   required to fit a fast changing function.  Too small a spread
%   means many neurons will be required to fit a smooth function,
%   and the network may not generalize well.  Call NEWRB with
%   different spreads to find the best value for a given problem.
%
%  Examples
%
%    Here we design a radial basis network given inputs P
%    and targets T.
%
%      P = [1 2 3];
%      T = [2.0 4.1 5.9];
%      net = srgtsnewrb(P,T);
%
%    Here the network is simulated for a new input.
%
%      P = 1.5;
%      Y = sim(net,P)
%
%  Algorithm
%
%    NEWRB creates a two layer network. The first layer has RADBAS
%    neurons, and calculates its weighted inputs with DIST, and
%    its net input with NETPROD.  The second layer has PURELIN neurons,
%    calculates its weighted input with DOTPROD and its net inputs with
%    NETSUM. Both layers have biases.
%
%    Initially the RADBAS layer has no neurons.  The following steps
%    are repeated until the network's mean squared error falls below GOAL
%   or the maximum number of neurons are reached:
%    1) The network is simulated
%    2) The input vector with the greatest error is found
%    3) A RADBAS neuron is added with weights equal to that vector.
%    4) The PURELIN layer weights are redesigned to minimize error.

display = 0;
% switch displayLevel
%     case { 'iter' , 'diagnose' }
%         display = 1;
% end

if nargin < 2
  net = newnet('newrb');
%   tr = [];
  return
end

% Defaults
if nargin < 3, goal = 0; end
if nargin < 4, spread = 1; end
if nargin < 6, df = 25; end

% Error checks
if (~isa(p,'double')) | (~isreal(p)) | (length(p) == 0)
  error('Inputs are not a non-empty real matrix.')
end
if (~isa(t,'double')) | (~isreal(t)) | (length(t) == 0)
  error('Targets are not a non-empty real matrix.')
end
if (size(p,2) ~= size(t,2))
  error('Inputs and Targets have different numbers of columns.')
end
if (~isa(goal,'double')) | ~isreal(goal) | any(size(goal) ~= 1) | (goal < 0)
  error('Performance goal is not a positive or zero real value.')
end
if (~isa(spread,'double')) | ~isreal(spread) | any(size(spread) ~= 1) | (spread < 0)
  error('Spread is not a positive or zero real value.')
end
if (~isa(df,'double')) | ~isreal(df) | any(size(df) ~= 1) | (df < 1) | (round(df) ~= df)
  error('Display frequency is not a positive integer.')
end

% More defaults
Q = size(p,2);
if nargin < 5, mn = Q; end

% More error checking
if (~isa(mn,'double')) | ~isreal(mn) | any(size(mn) ~= 1) | (mn < 1) | (round(mn) ~= mn)
  error('Maximum neurons is not a positive integer.')
end


% Dimensions
R = size(p,1);
S2 = size(t,1);

% Architecture
net = network(1,2,[1;1],[1; 0],[0 0;1 0],[0 1]);

% Simulation
net.inputs{1}.size = R;
net.layers{1}.size = 0;
net.inputWeights{1,1}.weightFcn = 'dist';
net.layers{1}.netInputFcn = 'netprod';
net.layers{1}.transferFcn = 'radbas';
net.layers{2}.size = S2;

% Performance
net.performFcn = 'mse';

% Design Weights and Bias Values
% [w1,b1,w2,b2,tr] = designrb(p,t,goal,spread,mn,df,display);
[w1,b1,w2,b2] = designrb(p,t,goal,spread,mn,df,display);

net.layers{1}.size = length(b1);
net.b{1} = b1;
net.iw{1,1} = w1;
net.b{2} = b2;
net.lw{2,1} = w2;

%======================================================
% function [w1,b1,w2,b2,tr] = designrb(p,t,eg,sp,mn,df,display)
function [w1,b1,w2,b2] = designrb(p,t,eg,sp,mn,df,display)

[r,q] = size(p);
[s2,q] = size(t);
b = sqrt(-log(.5))/sp;

% RADIAL BASIS LAYER OUTPUTS
P = radbas(dist(p',p)*b);
PP = sum(P.*P)';
d = t';
dd = sum(d.*d)';

% CALCULATE "ERRORS" ASSOCIATED WITH VECTORS
e = ((P' * d)' .^ 2) ./ (dd * PP');

% PICK VECTOR WITH MOST "ERROR"
pick = findLargeColumn(e);
used = [];
left = 1:q;
W = P(:,pick);
P(:,pick) = []; PP(pick,:) = [];
e(:,pick) = [];
used = [used left(pick)];
left(pick) = [];

% CALCULATE ACTUAL ERROR
w1 = p(:,used)';
a1 = radbas(dist(w1,p)*b);
[w2,b2] = solvelin2(a1,t);
a2 = w2*a1 + b2*ones(1,q);
sse = sumsqr(t-a2);

% % Start
% tr = newtr(mn,'perf');
% tr.perf(1) = sumsqr(t);
% tr.perf(2) = sse;

% if display
% 
%     if isfinite(df)
%         fprintf('NEWRB, neurons = 0, SSE = %g\n',sse);
%     end
% 
% end

% flag_stop = 0;

for k = 2:mn
  
  % CALCULATE "ERRORS" ASSOCIATED WITH VECTORS
  wj = W(:,k-1);
  a = wj' * P / (wj'*wj);
  P = P - wj * a;
  PP = sum(P.*P)';
  e = ((P' * d)' .^ 2) ./ (dd * PP');

  % PICK VECTOR WITH MOST "ERROR"
  pick = findLargeColumn(e);
  W = [W, P(:,pick)];
  P(:,pick) = []; PP(pick,:) = [];
  e(:,pick) = [];
  used = [used left(pick)];
  left(pick) = [];

  % CALCULATE ACTUAL ERROR
  w1 = p(:,used)';
  a1 = radbas(dist(w1,p)*b);
  [w2,b2] = solvelin2(a1,t);
  a2 = w2*a1 + b2*ones(1,q);
  sse = sumsqr(t-a2);
  
% %   % PROGRESS
% %   tr.perf(k+1) = sse;
% %   
% %   % DISPLAY
% %   if display
% %       if isfinite(df) && (~rem(k,df))
% %           fprintf('NEWRB, neurons = %g, SSE = %g\n',k,sse);
% %           flag_stop=plotperf(tr,eg,'NEWRB',k);
% %       end
% %   end
  % CHECK ERROR
  if (sse < eg), break, end
%   if (flag_stop), break, end

end

[S1,R] = size(w1);
b1 = ones(S1,1)*b;

% % Finish
% tr = cliptr(tr,k);

%======================================================
function i = findLargeColumn(m)

replace = find(isnan(m));
m(replace) = zeros(size(replace));

m = sum(m .^ 2,1);
i = find(m == max(m));
i = i(1);

%======================================================

function [w,b] = solvelin2(p,t)

if nargout <= 1
  w= t/p;
else
  [pr,pc] = size(p);
  x = t/[p; ones(1,pc)];
  w = x(:,1:pr);
  b = x(:,pr+1);
end

%======================================================
