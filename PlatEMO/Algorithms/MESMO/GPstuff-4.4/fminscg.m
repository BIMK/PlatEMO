function [x, fval, exitflag, output, grad] = fminscg(fun, x, opt)
%FMINSCG  Scaled conjugate gradient optimization
%
%  Description
%    X = FMINSCG(FUN, X0) starts at X0 and attempts to find a local
%    minimizer X of the function FUN. FUN accepts input X and
%    returns a scalar function value F and its scalar or vector
%    gradient G evaluated at X. X0 can be a scalar or vector
%
%    X = FMINSCG(FUN, X0, OPTIONS) minimizes with the default
%    optimization parameters replaced by values in the structure
%    OPTIONS, an argument created with the OPTIMSET function. See
%    OPTIMSET for details. Used options are Display, TolX, TolFun,
%    DerivativeCheck and MaxIter.
%
%    [X,FVAL] = FMINSCG(FUN,X0,...) returns the value of the objective 
%    function FUN at the solution X.
%
%    [X,FVAL,EXITFLAG] = FMINSCG(FUN,X0,...) returns an EXITFLAG that 
%    describes the exit condition of FMINSCG. Possible values of EXITFLAG 
%    and the corresponding exit conditions are listed below. See the
%    documentation for a complete description.
%      1  Magnitude of gradient small enough. 
%      2  Change in X too small.
%      3  Change in objective function too small.
%      0  Too many function evaluations or iterations.
%
%    [X,FVAL,EXITFLAG,OUTPUT] = FMINSCG(FUN,X0,...) returns a
%    structure OUTPUT with the number of iterations taken in
%    OUTPUT.iterations, the number of function evaluations in
%    OUTPUT.funcCount, the algorithm used in OUTPUT.algorithm,
%    function values for each iteration in OUTPUT.f and X for each
%    iteration in OUTPUT.x.
%
%    [X,FVAL,EXITFLAG,OUTPUT,GRAD] = FMINSCG(FUN,X0,...) returns the value 
%    of the gradient of FUN at the solution X.
%
%  Options
%    Display
%      - 'off' displays no output
%      - 'iter' displays output at each iteration, and gives the
%               default exit message.
%      - 'notify' displays output only if the function does not
%                 converge, and gives the default exit message.
%      - 'final' (default) displays just the final output, and
%                 gives the default exit message.
%    TolFun
%      Termination tolerance on the function value, a positive
%      scalar. The default is 1e-6.
%    TolX
%      Termination tolerance on x, a positive scalar. The default
%      value is 1e-6.
%    DerivativeCheck
%      Compare user-supplied derivatives (gradient of objective) to
%      finite-differencing derivatives. The choices are 'on' or the
%      default 'off'.
%    MaxIter
%      Maximum number of iterations allowed, a positive integer. 
%      The default value is 400.
%    lambda
%      Initial scale parameter for optimizer. The default value is 10
%    lambdalim
%      Limit for the scale paramer. The defualt value is 1e101.
%
%  See also OPTIMSET
%
% Copyright (c) 1996,1997 Christopher M Bishop, Ian T Nabney
% Copyright (c) 2005,2010,2012 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

tic
% Set empty options to default values
defaultopt = struct( ...
    'DerivativeCheck','off', ...   
    'Display','final', ...
    'MaxIter',400, ...
    'TolFun',1e-6, ...
    'TolX',1e-6, ...
    'lambda', 10, ...
    'lambdalim', 1e20, ...
    'GradConstr', 'off'); 

% If just 'defaults' passed in, return the default options in X
if nargin==1 && nargout <= 1 && isequal(fun,'defaults')
   x = defaultopt;
   return
end

if nargin < 3, opt=[]; end 

switch optimget(opt,'Display',defaultopt,'fast');
  case 'off'
    display=0;
  case 'notify'
    display=1;
  case 'final'
    display=2;
  case 'iter'
    display=3;
  otherwise
    display=2;
end
maxiter = optimget(opt,'MaxIter',defaultopt,'fast');
tolfun = optimget(opt,'TolFun',defaultopt,'fast');
tolx = optimget(opt,'TolX',defaultopt,'fast');
lambda = optimget(opt,'lambda', defaultopt, 'fast');
lambdalim = optimget(opt, 'lambdalim', defaultopt, 'fast');
GradConstr = optimget(opt, 'GradConstr', defaultopt, 'fast');

nparams = length(x);

%  Check gradients
if isequal(optimget(opt,'DerivativeCheck',defaultopt,'fast'),'on');
  derivativecheck(x, fun);
end

if (display >= 3)
  if isequal(GradConstr,'on')
    fprintf('  Iteration  Func-count  Grad-count    f(x)    Lambda\n');
  else
    fprintf('  Iteration  Func-count      f(x)      Lambda\n');
  end
end

sigma0 = 1.0e-4;
[fold,gradold] = fun(x); % Initial function value and gradient
funcCount=1;
gradCount=1;
gradnew = gradold;
d = - gradnew;           % Initial search direction.
success = 1;             % Force calculation of directional derivs.
nsuccess = 0;            % nsuccess counts number of successes.
lambdamin = 1.0e-15; 
lambdamax = 1.0e100;

if (display >= 3)
  if isequal(GradConstr,'on')
    fprintf('  %5.0f      %5.0f       %5.0f  %12.5g    \n',0,funcCount,gradCount,fold);
  else
    fprintf('  %5.0f      %5.0f    %12.5g    \n',0,funcCount,fold);
  end
end

j = 1;                   % j counts number of iterations.
if nargout >= 4
  output.f(j, :) = fold;
  output.x(j, :) = x;
  output.algorithm='fminscg';
end

% Main optimization loop.
while (j <= maxiter)

  % Calculate first and second directional derivatives.
  if (success == 1)
    mu = d*gradnew';
    if (mu >= 0)
      d = - gradnew;
      mu = d*gradnew';
    end
    kappa = d*d';
    if kappa < eps
      if (display >= 2)
        disp('Gradient smaller than eps');
      end
      exitflag=1;
      if nargin>4
        output.funcCount=funcCount;
      end
      return
    end
    sigma = sigma0/sqrt(kappa);
    xplus = x + sigma*d;
    [tmp,gplus] = fun(xplus);
    funcCount=funcCount+1;
    gradCount=gradCount+1;
    gamma = (d*(gplus' - gradnew'))/sigma;
  end

  % Increase effective curvature and evaluate step size alpha.
  delta = gamma + lambda*kappa;
  if (delta <= 0)
    delta = lambda*kappa;
    lambda = lambda - gamma/kappa;
  end
  alpha = - mu/delta;
  
  % Calculate the comparison ratio.
  xnew = x + alpha*d;
  if isequal(GradConstr,'on')
    fnew = fun(xnew);
    funcCount=funcCount+1;
  else
    [fnew,gnew] = fun(xnew);
    funcCount=funcCount+1;
    gradCount=gradCount+1;
  end
  while isinf(fnew) || isnan(fnew)
    warning('Function value at xnew not finite or a number')
    lambda = min(4.0*lambda, lambdamax);
    delta = gamma + lambda*kappa;
    if (delta <= 0)
      delta = lambda*kappa;
      lambda = lambda - gamma/kappa;
    end
    alpha = - mu/delta;
    xnew = x + alpha*d;
    if isequal(GradConstr,'on')
      fnew = fun(xnew);
      funcCount=funcCount+1;
    else
      [fnew,gnew] = fun(xnew);
      funcCount=funcCount+1;
      gradCount=gradCount+1;
    end
  end
  Delta = 2*(fnew - fold)/(alpha*mu);
  if (Delta  >= 0)
    success = 1;
    nsuccess = nsuccess + 1;
    x = xnew;
    fnow = fnew;
  else
    success = 0;
    fnow = fold;
  end

  if nargout >= 4
    % Store relevant variables
    output.f(j) = fnow;               % Current function value
    output.x(j,:) = x;      % Current position
  end    
  if display >= 3
    if rem(j,20)==0
      if isequal(GradConstr,'on')
        fprintf('  Iteration  Func-count  Grad-count    f(x)    Lambda\n');
      else
        fprintf('  Iteration  Func-count      f(x)    Lambda\n');
      end
    end
    if isequal(GradConstr,'on')
      fprintf('  %5.0f      %5.0f       %5.0f  %12.6g   %6.3g\n',j,funcCount,gradCount,fnow,lambda);
    else
      fprintf('  %5.0f      %5.0f    %12.6g   %6.3g\n',j,funcCount,fnow,lambda);
    end
  end
  
  if (success == 1)
    
    % Test for termination
    if max(abs(alpha*d)) < tolx
      if nargin <5
        fval=fnew;
      else
        [fval,grad]=fun(x);
        funcCount=funcCount+1;
        gradCount=gradCount+1;
      end
      exitflag=2;
      if nargin>4
        output.funcCount=funcCount;
      end
      if (display >= 2)
        if isequal(GradConstr,'on')
          fprintf(' TolX reached. Func-count %d. Grad-count %d. Final f(x)=%6.6g. Elapsed time %.2f\n',funcCount,gradCount,fval,toc)
        else
          fprintf(' TolX reached. Func-count %d. Final f(x)=%6.6g. Elapsed time %.2f\n',funcCount,fval,toc)
        end
      end
      return
      
    elseif max(abs(fnew-fold)) < tolfun
      if nargin <5
        fval=fnew;
      else
        [fval,grad]=fun(x);
        funcCount=funcCount+1;
        gradCount=gradCount+1;
      end
      exitflag=3;
      if nargin>4
        output.funcCount=funcCount;
      end
      if (display >= 2)
        if isequal(GradConstr,'on')
          fprintf(' TolFun reached. Func-count %d. Grad-count %d. Final f(x)=%6.6g. Elapsed time %.2f\n',funcCount,gradCount,fval,toc)
        else
          fprintf(' TolFun reached. Func-count %d. Final f(x)=%6.6g. Elapsed time %.2f\n',funcCount,fval,toc)
        end
      end
      return

    else
      % Update variables for new position
      if isequal(GradConstr,'on')
        fold = fnew;
        gradold = gradnew;
        [fval,gradnew] = feval(fun, x);
        funcCount=funcCount+1;
        gradCount=gradCount+1;
      else
        fold = fnew;
        fval = fnew;
        gradold = gradnew;
        gradnew = gnew;
      end
      % If the gradient is zero then we are done.
      if (gradnew*gradnew' < eps) && all(isreal(gradnew))
        grad=gradnew;
        exitflag=1;
        if nargin>4
          output.funcCount=funcCount;
        end
      if (display >= 2)
        if isequal(GradConstr,'on')
          fprintf(' Gradient smaller than eps. Func-count %d. Grad-count %d. Final f(x)=%6.6g. Elapsed time %.2f\n',funcCount,gradCount,fval,toc)
        else
          fprintf(' Gradient smaller than eps. Func-count %d. Final f(x)=%6.6g. Elapsed time %.2f\n',funcCount,fval,toc)
        end
      end
        return
      end
    end
  end

  % Adjust lambda according to comparison ratio.
  if (Delta < 0.25)
    lambda = min(4.0*lambda, lambdamax);
  end
  if (Delta > 0.75)
    lambda = max(0.5*lambda, lambdamin);
  end
  
  % If scale parameter is at its limit, stop optimization
  if lambda >= lambdalim
    exitflag=0;
    if (display >= 1)
      disp(['Warning: Optimization stopped because lambda parameter reached limit';...
            '         Check that the analytic gradients are correct!             ']);
    end
    grad=gradnew;
    fval=fnew;
    return
  end

  % Update search direction using Polak-Ribiere formula, or re-start 
  % in direction of negative gradient after nparams steps.
  if (nsuccess == nparams)
    d = -gradnew;
    nsuccess = 0;
  else
    if (success == 1)
      beta = (gradold - gradnew)*gradnew'/(mu);
      d = beta*d - gradnew;
    end
  end
  
  j = j + 1;
  output.iterations=j;
  
end

% If we get here, then we haven't terminated in the given number of 
% iterations.
exitflag=0;
if (display >= 1)
  disp('Maximum number of iterations has been exceeded.');
end
grad=gradnew;
fval=fnew;
if (display >= 2)
  if isequal(GradConstr,'on')
    fprintf(' Func-count %d. Grad-count %d. Final f(x)=%6.6g. Elapsed time %.2f\n',funcCount,gradCount,fval,toc)
  else
    fprintf(' Func-count %d. Final f(x)=%6.6g. Elapsed time %.2f\n',funcCount,fval,toc)
  end
end
