function obj=Infill_Standard_EI(x, Kriging_model, f_min)
%--------------------------------------------------------------------------
% the DACE toolbox of  Lophaven et al. (2002)  is used to predict value
%--------------------------------------------------------------------------
% REFERENCES
%Lophaven SN, Nielsen HB, and Sodergaard J, DACE - A MATLAB Kriging
%Toolbox, Technical Report IMM-TR-2002-12, Informatics and Mathematical
%Modelling, Technical University of Denmark, 2002.
%Available at: http://www2.imm.dtu.dk/~hbn/dace/.
%--------------------------------------------------------------------------
% get the Kriging prediction and variance
[y,mse] = predictor(x,Kriging_model);
s=sqrt(max(0,mse));
% calcuate the EI value
EI=(f_min-y).*Gaussian_CDF((f_min-y)./s)+s.*Gaussian_PDF((f_min-y)./s);

% the genetic algorithm tries to minimize the objective
obj=-EI;

end

function res=Gaussian_PDF(x)
res=1/sqrt(2*pi)*exp(-x.^2/2);
end

function x=Gaussian_CDF(x)
  x=0.5*(1+erf(x/sqrt(2)));
end



