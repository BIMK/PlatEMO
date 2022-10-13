% KEVALSRESET - Resets the kernel evaluations counter 
%
% Syntax: kevalsreset();
%
% Version 3.22e -- Comments to diehl@alumni.cmu.edu
%

function kevalsreset();

global kernel_evals;

kernel_evals = 0;
