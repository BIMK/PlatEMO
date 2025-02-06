function y_i = SolutionSelection(Problem, Population, P, model, N_s, i)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Haoran Gu

    Prnd = nchoosek(P,2);
    Prnd = [Prnd;Prnd];
    Prnd(N_s/2+1:N_s,[1,2]) = Prnd(N_s/2+1:N_s,[2,1]);
    for r = 1 : N_s
        candidate(r,:) = OperatorDE(Problem, Population(i).dec, Population(Prnd(r,1)).dec, Population(Prnd(r,2)).dec);
    end
    c       = model.PredictClass(candidate);
    [~,pos] = min(c);
    y_i     = candidate(pos,:);
end