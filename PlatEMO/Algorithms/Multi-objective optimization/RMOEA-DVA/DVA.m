function  [HR,LR] = DVA(Problem,Population,nDVA,theta)
% Decision variable assortment

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [N,D] = size(Population.decs);
    sd    = zeros(nDVA,D);
    for i = 1 : D
        Decs = Population(randperm(N,nDVA)).decs;
        for j = 1 : nDVA
            % Perturb the i-th variable of the j-th solution
            delta     = Problem.delta*(Problem.upper(i)-Problem.lower(i));
            PDec      = repmat(Decs(j,:),Problem.H,1);
            PDec(:,i) = PDec(:,i) + 2*delta*rand(Problem.H,1) - delta;
            NewP      = Problem.Evaluation(PDec);
            FrontNo   = NDSort(NewP.objs,inf);
            % Calculate the variance
            sd(j,i) = std(FrontNo);
        end
    end
    sd = (sd-min(sd(:)))./(max(sd(:))-min(sd(:)));
    M  = mean(sd,1);
    HR = find(M>theta);
    LR = find(M<=theta);
end