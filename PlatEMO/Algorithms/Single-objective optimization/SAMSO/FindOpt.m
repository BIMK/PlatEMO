function  opt = FindOpt(model,Population,BU,BD)
% Find the minimum of the surrogate

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [~,I]  = sort(Population.objs,'ascend');
    init   = 1;
    preObj = rbf_predict(model,Population.decs,Population(I(init)).decs);
    while isnan(preObj)
        init   = init + 1;
        preObj = rbf_predict(model,Population.decs,Population(I(init)).decs);
    end
    opt = fmincon(@(x)rbf_predict(model,Population.decs,x),Population(I(init)).decs,[],[],[],[],BD,BU,[],optimoptions('fmincon','Display','off'));
end