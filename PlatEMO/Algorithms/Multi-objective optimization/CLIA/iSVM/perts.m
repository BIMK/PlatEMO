% PERTS - Returns the number of perturbations 
%
% Syntax: perturbations = perts();
%
% Version 3.22e -- Comments to diehl@alumni.cmu.edu
%

function perturbations = perts();

global perturbations;

if (isempty(perturbations))
   perturbations = 0;
end;
