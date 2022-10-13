% KEVALS - Returns the number of kernel evaluations 
%
% Syntax: kernel_evals = kevals();
%
% Version 3.22e -- Comments to diehl@alumni.cmu.edu
%

function kernel_evals = kevals();

global kernel_evals;

if (isempty(kernel_evals))
   kernel_evals = 0;
end;
