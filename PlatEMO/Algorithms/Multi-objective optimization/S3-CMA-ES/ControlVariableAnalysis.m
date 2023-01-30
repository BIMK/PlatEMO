function [PV,DV] = ControlVariableAnalysis(Problem,NCA)
% Control variable analysis

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Huangke Chen
    
    Fno = zeros(1,Problem.D);
    for i = 1 : Problem.D
        x       = 0.2*ones(1,Problem.D).*(Problem.upper-Problem.lower) + Problem.lower;
        S       = repmat(x,NCA,1);
        inter   = (0.95-0.05)/(NCA-1);
        tempA   = 0.05:inter:0.96;
        S(:, i) = tempA'*(Problem.upper(i)-Problem.lower(i)) + Problem.lower(i);
        S          = Problem.Evaluation(S);
        [~,MaxFNo] = NDSort(S.objs,inf);
        Fno(i) = MaxFNo;      
    end
    [~,I] = sort(Fno);
    PV    = sort(I(1:Problem.M-1));
    DV    = sort(I(Problem.M:end));
end