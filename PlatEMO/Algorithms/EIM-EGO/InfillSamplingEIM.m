function PopDec = InfillSamplingEIM(Global,KrigingModel,PopObjScaled,InfillCriterionIndex)
% Solution update in MOEGO, where a solution with the best expected
% improvement matrix value is selected

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Dawei Zhan

    % the number of design variables
    D = Global.D;
    % number of objective functions
    M = Global.M;
    % the non-dominated front
    index = NDSort(PopObjScaled,1)==1;
    f = PopObjScaled(index',:);
    % number of non-dominated front points
    p = size(f,1);
    % population size of GA
    GAPopulationSize = 100;
    GAGeneration = 100;
    EIM_max    = inf;
    % the first GA generation, randomly generated
    Offspring   = repmat(Global.upper-Global.lower,GAPopulationSize,1).*lhsamp(GAPopulationSize,D)+repmat(Global.lower,GAPopulationSize,1);
    % the GA process for optimizing the EIM function
    for gen = 1 :  GAGeneration
        % the kriging prediction and varince
        u = zeros(GAPopulationSize,M);
        mse = zeros(GAPopulationSize,M);
        for i = 1 : M
            [u(:, i),mse(:, i)] = predictor(Offspring, KrigingModel{i});
        end
        s = sqrt(max(0,mse));
        % the EI matrix (three dimensional matrix)
        f_matrix  =  f .* ones(1,1,GAPopulationSize);
        u_matrix = reshape(u', 1, M, GAPopulationSize).* ones(p,1);
        s_matrix =  reshape(s', 1, M, GAPopulationSize).* ones(p,1);
        EI_matrix=(f_matrix-u_matrix).*Gaussian_CDF((f_matrix-u_matrix)./s_matrix)+s_matrix.*Gaussian_PDF((f_matrix-u_matrix)./s_matrix);
        switch InfillCriterionIndex
            case 1
                %  the Euclidean distance-based EI matrix criterion
                EIM = -reshape(min(sqrt(sum(EI_matrix.^2,2)), [], 1), GAPopulationSize, 1, 1);
            case 2
                %  the Maximin distance-based EI matrix criterion
                EIM = -reshape(min(max(EI_matrix, [], 2), [], 1), GAPopulationSize, 1, 1);
            case 3
                %  the Hypervolume-based EI matrix criterion
                EIM = -reshape(min(prod(1.1*ones(1,M)-f+EI_matrix,2)-prod(1.1*ones(1,M)-f,2), [], 1), GAPopulationSize, 1, 1);
        end
        [~,index] = sort(EIM,'ascend');
        if EIM(index(1)) < EIM_max
            Best = Offspring(index(1),:);
            EIM_max   = EIM(index(1));
        end
        Parent    = Offspring(index(1:ceil(GAPopulationSize/2)),:);
        Offspring = [GA(Parent(TournamentSelection(2,size(Parent,1),EIM(index(1:ceil(GAPopulationSize/2)))),:));
                     GA(Parent,{0,0,1,20})];
    end

    PopDec = Best;

end

function y=Gaussian_PDF(x)
    y=1/sqrt(2*pi)*exp(-x.^2/2);
end

function y=Gaussian_CDF(x)
	y=0.5*(1+erf(x/sqrt(2)));
end