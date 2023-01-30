function newDec = GenerateNew(Problem,Population,Lammda,Mu,CR)
% Generate offspring in SADE-Sammon

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    % Get population information
    Decs = Population.decs;
    Objs = Population.objs;
    
    % Sort current evaluated solution
    [~,index] = sort(Objs,'ascend');
    bestDec   = Decs(index(1),:);
    parentIdx2 = randperm(Lammda);
    parentIdx3 = randperm(Lammda);
    % Check parents
    while any(parentIdx2==parentIdx3)
        same = parentIdx2==parentIdx3;
        parentIdx3(same) = randperm(Lammda,sum(same));
    end
    
    % DE/best/1
    Parent1 = repmat(bestDec,Lammda,1);
    Parent2 = Decs(parentIdx2,:);
    Parent3 = Decs(parentIdx3,:);
    newDec  = OperatorDE(Problem,Parent1,Parent2,Parent3,{CR,Mu,1,20});
end