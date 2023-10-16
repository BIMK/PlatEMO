function outcome = Normalization(Population,Average,G)
% Normalization

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    fmax = max(Population.objs,[],1);
    fmin = min(Population.objs,[],1);
    Average(G,:) = mean((Population.objs-fmin)./(fmax-fmin));
    outcome = Average;
end