function TP = UpdateTP(Population,W)
% Update the target points

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    TP = [];
    [Subpopulation,Z] = Classification(Population,W);
    for i = 1 : length(Subpopulation)
        ASF = max((Subpopulation{i}.objs-repmat(Z,length(Subpopulation{i}),1))./repmat(W(i,:),length(Subpopulation{i}),1),[],2);
        [~,extreme] = min(ASF);
        TP = [TP,Subpopulation{i}(extreme)];
    end
end