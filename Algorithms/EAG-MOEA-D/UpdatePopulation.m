function Population = UpdatePopulation(Population,Offspring,offspringLoc,W,B)
% Update the population by MOEA/D

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    for i = 1 : length(Offspring)
        % Update the parents by weight sum approach
        g_old = sum(Population(B(offspringLoc(i),:)).objs.*W(B(offspringLoc(i),:),:),2);
        g_new = W(B(offspringLoc(i),:),:)*Offspring(i).obj';
        Population(B(offspringLoc(i),g_old>=g_new)) = Offspring(i);
    end
end