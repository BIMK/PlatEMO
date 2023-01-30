function Population = DecompositionSelection(Global,Population,associate,Cosinemax)
% The decomposition-based method environmental selection

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Shufen Qin
% E-mail: shufen.qin@stu.tyust.edu.cn

    np = length(Population);
    %% Normalization
    Obj = Population.objs;
    Obj = (Obj-repmat(min(Obj),np,1))./(repmat(max(Obj),np,1)-repmat(min(Obj),np,1));

    %% Select one solution for each reference vector
    list = unique(associate)';
    Next = zeros(length(list),1);
    t = 1;
    for i = list
        current = find(associate == i);
        dist = pdist2(Obj(current,:),zeros(1,Global.M),'Euclidean');
        Fan = Cosinemax(current)./dist;
        [~,best] = max(Fan);
        Next(t)  = current(best);
        t = t +1;
    end
    % Population for next generation
    Population = Population(Next);
end