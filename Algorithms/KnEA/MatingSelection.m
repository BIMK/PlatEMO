function MatingPool = MatingSelection(PopObj,FrontNo,KneePoints)
% The mating selection of KnEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Calculate the weighted distance of each solution
    Dis   = pdist2(PopObj,PopObj);
    Dis(logical(eye(length(Dis)))) = inf;
    Dis   = sort(Dis,2);
	Crowd = sum(Dis(1:3,:).*repmat((3:-1:1)',1,size(PopObj,1)));

    %% Binary tournament selection
    MatingPool = TournamentSelection(2,size(PopObj,1),FrontNo,-KneePoints,-Crowd);
end