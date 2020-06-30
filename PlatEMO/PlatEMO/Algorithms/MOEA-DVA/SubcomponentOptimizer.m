function Population = SubcomponentOptimizer(Population,Neighbour,indices)
% Subcomponent optimizer

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    for i = 1 : length(Population)
        if rand < 0.9
            P = Neighbour(i,randperm(size(Neighbour,2),2));
        else
            P = randperm(length(Population),2);
        end
        OffDec          = Population(i).dec;
        NewDec          = DE(OffDec,Population(P(1)).dec,Population(P(2)).dec,{1,0.5,length(OffDec)/length(indices)/2,20});
        OffDec(indices) = NewDec(indices);
        Offspring       = INDIVIDUAL(OffDec);
        if sum(Offspring.obj) < sum(Population(i).obj)
            Population(i) = Offspring;
        end
    end
end