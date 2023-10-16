function  [HR,LR] = DVC(Problem,Population,SN,PN,TN,theta)
% Decision variable classification

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

	[N,D]  = size(Population.decs);
    AllVal = zeros(TN,D);
    for T = 1 : TN
        Var = zeros(SN,D);
        for i = 1 : D
            a    = randperm(N,SN);
            Decs = Population(a).decs;
            Objs = Population(a).objs;
            for j = 1 : SN
                % Perturb the i-th variable of the j-th solution for PN times
                delta     = Problem.delta*(Problem.upper(i)-Problem.lower(i));
                PDec      = repmat(Decs(j,:),PN,1);
                PDec(:,i) = PDec(:,i) + 2*delta*rand(PN,1) - delta;
                NewP      = Problem.Evaluation(PDec);
                % Calculate the variance
                vcv      = sum(abs(NewP.objs-repmat(Objs(j,:),PN,1)),2);
                Var(j,i) = std(vcv,0,1)^2;
            end
        end
        AllVal(T,:) = mean(Var,1);
    end
    TVal = zeros(1,D);
    for i = 1 : D
        TVal(i) = numel(find(AllVal(:,i)<theta));
    end
    HR = find(TVal<=mean(TVal));	% Highly robustness-related variables
    LR = find(TVal>mean(TVal));     % Weakly robustness-related variables
end