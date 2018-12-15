function a = updateCovariance(a,xstep)
% Update the covariance of CMA model

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    % Constant of learning rate for evolution path
    cc = 2/(length(a.x)+2);
    % Constant of learning rate for covariance matrix
    ccov = 2/(length(a.x)^2+6);
    if a.psucc < 0.44
        % Update the evolution path
        a.pc = (1-cc)*a.pc + sqrt(cc*(2-cc))*xstep;
        % Update the covariance matrix
        a.C = (1-ccov)*a.C + ccov*a.pc'*a.pc;
    else
        % Update the evolution path
        a.pc = (1-cc)*a.pc;
        % Update the covariance matrix
        a.C = (1-ccov)*a.C + ccov*(a.pc'*a.pc+cc*(2-cc)*a.C);
    end
end