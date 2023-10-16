function Population = EnvironmentalSelection(Population,Offspring)
% Environmental selection of FROFI

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Population update
    replace1 = FitnessSingle(Population) > FitnessSingle(Offspring);
    replace2 = ~replace1 & Population.objs > Offspring.objs;
    Population(replace1) = Offspring(replace1);
    Archive = Offspring(replace2);
    
    %% Replacement mechanism
    Nf = round(length(Population)/round(max(5,length(Population(1).dec)/2)));
    [~,rank]   = sort(Population.objs,'descend');
    Population = Population(rank);
    for i = 1 : floor(length(Population)/Nf)
        if isempty(Archive)
            break;
        else
            current   = (i-1)*Nf+1 : i*Nf;
            [~,worst] = max(sum(max(0,Population(current).cons),2));
            [~,best]  = min(sum(max(0,Archive.cons),2));
            if Archive(best).obj < Population(current(worst)).obj
                Population(current(worst)) = Archive(best);
                Archive(best) = [];
            end
        end
    end
end