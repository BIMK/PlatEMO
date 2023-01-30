function ProposedPoint = GAOptimize(model,PopSize,Decs,BU,BD)
% GA optimization in SACCEAMII

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    D = size(Decs,2);
    PopDec   = rand(PopSize,D).*repmat(BU-BD,PopSize,1) + repmat(BD,PopSize,1);
    PopObj   = rbf_predict(model,Decs,PopDec);
    Population = [PopDec,PopObj];
    maxIter  = 100;
    iter = 0;
    while iter < maxIter
        MatingPool = TournamentSelection(2,PopSize,Population(:,D+1));
        OffDec     = GAOperator(Population(MatingPool,1:D),BU,BD);
        OffObj     = rbf_predict(model,Decs,OffDec);
        Population = [Population;OffDec,OffObj];
        [~,rank]   = sort(Population(:,D+1),'ascend');
        Population = Population(rank(1:PopSize),:);
        iter = iter + 1;
    end
    [~,best] = min(Population(:,D+1));
    ProposedPoint = Population(best,1:D);
end