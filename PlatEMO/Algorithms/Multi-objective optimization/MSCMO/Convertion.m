function outcome = Convertion(Average,G,g,last_gen,change_threshold)
% Conversion condition for enter next stage

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    if (G-g > last_gen) && (max(abs(Average(G,:)-Average(G-last_gen,:)))<=change_threshold)
        outcome = true;
    else
        outcome = false;
    end
end