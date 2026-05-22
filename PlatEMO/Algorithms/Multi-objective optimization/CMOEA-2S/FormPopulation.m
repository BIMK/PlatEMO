function [Population, PopulationD, PopIndex] = FormPopulation(Archive, Action, PopSize)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    ArcSize = numel(Archive);

    %% Form population(s)
    if ArcSize < PopSize
        PopIndex    = 1 : 1 : ArcSize;
        Population  = Archive;
        PopulationD = Archive;
    else
        if Action == 1
            [Population,~, PopIndex] = AdaSelection(Archive,PopSize,true);
            PopulationD = [];
        elseif Action == 2
            if rand()<0.1
                [~,~, PopIndex] = AdaSelection(Archive,PopSize,true);
                PopIndex        = ~PopIndex;
                PopulationD     = Archive(PopIndex);
            else
                [PopulationD,~, PopIndex] = AdaSelection(Archive,PopSize,false);
            end
            Population = [];
        end
    end
end