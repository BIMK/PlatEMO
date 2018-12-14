function Next = Classification(Global,Population,W,Z)
% Classify solutions into sub-regions

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Calculate the cosine value of each solution to each vector
    % Note that here the value of each solution on each vector is not the
    % cosine value which is proposed in the paper
    PopObj = Population.objs - repmat(Z,length(Population),1);
    Value  = PopObj*W';
    [~,P]  = max(Value,[],2);
    
    %% Select one solution for each sub-region
    Next = Population(1:size(W,1));
    for i = 1 : size(W,1)
        Current = find(P==i);
        if isempty(Current)
            Next(i) = Global.Initialization(1);
        else
            ND       = find(NDSort(Population(Current).objs,1)==1);
            [~,best] = max(PopObj(Current(ND),:)*W(i,:)'./sum(PopObj(Current(ND),:).^2,2).^0.6);
            Next(i)  = Population(Current(ND(best)));
        end
    end
end