function [Arch,score] = Initialization(Population,Dec,Mask,FrontNo)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [Next]     =  find(FrontNo==1);
    Population = Population(Next).objs;
    Dec        = Dec(Next,:);
    Mask       = Mask(Next,:);
    score      = sum(Mask,1)/size(Mask,1);
    Arch       = archives(Population,Dec,Mask);
end