function Offspring = Operator(Population,K)
% Generate offsprings by the models
% K --- 5 --- Number of clusters

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is modified from the code in
% http://dces.essex.ac.uk/staff/zhang/IntrotoResearch/RegEDA.htm

    %% Parameter setting
    if nargin < 2
        K = 5;
    end
    PopDec = Population.decs;
    [N,D]  = size(PopDec);
    M      = length(Population(1).obj);
    
    %% Modeling
    [Model,probability] = LocalPCA(PopDec,M,K);

    %% Reproduction
    OffspringDec = zeros(N,D);
    % Generate new trial solutions one by one
    for i = 1 : N
        % Select one cluster by Roulette-wheel selection
        k = find(rand<=probability,1);
        % Generate one offspring
        if ~isempty(Model(k).eVector)
            lower = Model(k).a - 0.25*(Model(k).b-Model(k).a);
            upper = Model(k).b + 0.25*(Model(k).b-Model(k).a);
            trial = rand(1,M-1).*(upper-lower) + lower;
            sigma = sum(abs(Model(k).eValue(M:D)))/(D-M+1);
            OffspringDec(i,:) = Model(k).mean + trial*Model(k).eVector(:,1:M-1)' + randn(1,D)*sqrt(sigma);
        else
            OffspringDec(i,:) = Model(k).mean + randn(1,D);
        end
    end
    Offspring = INDIVIDUAL(OffspringDec);
end