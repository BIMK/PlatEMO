% PLOT2DKM - For a 2-D binary classification problem, plot2dkm plots the data, 
%            the margin and error vectors and contours of constant margin
%            for the SVM classifier in memory.
%
% Syntax: plot2dkm
%
% Version 3.22e -- Comments to diehl@alumni.cmu.edu
%

function plot2dkm

% flags for example state
MARGIN    = 1;
ERROR     = 2;
RESERVE   = 3;
UNLEARNED = 4;

% define global variables 
global ind;   % cell array containing indices of margin, error, reserve and unlearned vectors
global X;     % matrix of margin, error, reserve and unlearned vectors stored columnwise
global y;     % column vector of class labels (-1/+1) for margin, error, reserve and unlearned vectors

% plot examples with label -1
figure;
indn1 = find(y == -1);
scatter(X(1,indn1),X(2,indn1),40,'b','filled');
hold on;

% plot examples with label +1
ind1 = find(y == 1);
scatter(X(1,ind1),X(2,ind1),40,'r');

% plot margin vectors
scatter(X(1,ind{MARGIN}),X(2,ind{MARGIN}),120,'k');
scatter(X(1,ind{MARGIN}),X(2,ind{MARGIN}),150,'k');
scatter(X(1,ind{MARGIN}),X(2,ind{MARGIN}),200,'k');

% plot error vectors
scatter(X(1,ind{ERROR}),X(2,ind{ERROR}),120,'k');
scatter(X(1,ind{ERROR}),X(2,ind{ERROR}),150,'k');
scatter(X(1,ind{ERROR}),X(2,ind{ERROR}),200,'k');
scatter(X(1,ind{ERROR}),X(2,ind{ERROR}),120,'k','x');
scatter(X(1,ind{ERROR}),X(2,ind{ERROR}),150,'k','x');
scatter(X(1,ind{ERROR}),X(2,ind{ERROR}),200,'k','x');

% draw margin band
xl = xlim;
yl = ylim;
pd = min(xl(2)-xl(1),yl(2)-yl(1))/100;
x_range = xl(1):pd:xl(2);
y_range = yl(1):pd:yl(2);
f = zeros(length(y_range),length(x_range));
i = 1;
for xp = x_range
   j = 1;
   for yp = y_range
      f(j,i) = svmeval([xp ; yp]);
      j = j + 1;
   end;
   i = i + 1;
end;
contour(x_range,y_range,f,[-1 1]);





