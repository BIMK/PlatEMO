function TP = UpdateTP(Population,W)
% Update the target points

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
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