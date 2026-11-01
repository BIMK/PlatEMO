function [Population,Dec,Mask,Fitness,Importance] = MASSRMARLACSInitialize(Problem)
%MASSRMARLACSInitialize - Estimate variable importance and initialize search.
%
%   Single-variable probes estimate importance. Real-valued problems use
%   five probe sets and binary problems use one. Retain the first probe set
%   with the initial masked population for environmental selection.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    D          = Problem.D;
    Importance = zeros(1,D);
    Temp       = [];
    TDec       = [];
    TMask      = [];
    if any(Problem.encoding~=4)
        repeats = 5;
    else
        repeats = 1;
    end

    % Estimate the importance of activating each variable
    for sample = 1 : repeats
        value                      = Problem.lower + rand(1,D).*(Problem.upper-Problem.lower);
        value(Problem.encoding==4) = 1;
        probes                     = repmat(value,D,1);
        active                     = eye(D) > 0;
        P                          = Problem.Evaluation(probes.*active);
        Importance                 = Importance + sum(P.objs,2)';
        if sample == 1
            Temp  = P;
            TDec  = probes;
            TMask = active;
        end
    end

    % Create latent decisions and tournament-guided masks
    Dec                        = Problem.lower + rand(Problem.N,D).*(Problem.upper-Problem.lower);
    Dec(:,Problem.encoding==4) = 1;
    Mask                       = false(Problem.N,D);
    for i = 1 : Problem.N
        Mask(i,TournamentSelection(2,max(1,ceil(rand*D)),Importance)) = true;
    end
    P                             = Problem.Evaluation(Dec.*Mask);
    [Population,Dec,Mask,Fitness] = MASSRMARLACSSelection([Temp,P],[TDec;Dec],[TMask;Mask],Problem.N);
end