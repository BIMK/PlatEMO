function Population = UpdatePopulation(Population,Offspring)
% Update the population

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N = length(Population);

    %% Insert the offspring into the population by epsilon-dominance
    if ~any(all(Population.objs<=repmat(Offspring.obj,N,1),2))
        Dominate = find(all(repmat(Offspring.obj,N,1)<=Population.objs,2));
        if ~isempty(Dominate)
            k = randi(length(Dominate));
            Population(Dominate(k)) = Offspring;
        else
            k = randi(N);
            Population(k) = Offspring;
        end
    end
end