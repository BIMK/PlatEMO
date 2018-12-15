function [wz,AA,RA] = Interaction(PopObj,Point)
% Identify the preferred point and all the others are treated as
% non-preferred points

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Identify the unique and non-dominated solutions
    PopObj = unique(PopObj,'rows');
    PopObj = PopObj(NDSort(PopObj,1)==1,:);
    
    %% Identify the preferred solution
    ideal = min(PopObj,[],1);
    % Calculate the Tchebycheff function value of each solution on Point
    Fitness = max((PopObj-repmat(ideal,size(PopObj,1),1))./repmat(Point,size(PopObj,1),1),[],2);
    % The one having the minimal function value is treated as the preferred
    % point
    [~,prefer] = min(Fitness);
    AA = PopObj(prefer,:);
    RA = PopObj([1:prefer-1,prefer+1:end],:);
    
    %% Calculate the weight distribution function value
    RefPoint = max(PopObj,[],1) + 0.1;
    wz = [0,1,1+CalWHV(AA,RefPoint,ones(1,size(AA,1)))/CalWHV(RA,RefPoint,ones(1,size(RA,1)))];
end