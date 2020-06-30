function MatingPool = MatingSelection(PopObj,Rank)
% The mating selection of one-by-one EA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    %% Calculate the density of each solution
    d  = pdist2(PopObj,PopObj,'cosine');
    d  = sort(d,2);
    dk = 1./(sum(d(:,2:ceil(end/10)),2)+1);

    %% Binary tournament selection
    MatingPool = TournamentSelection(2,size(PopObj,1),Rank,dk);
end