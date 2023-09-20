function index = KrigingSelection(PopObj,V)
% The environmental selection operator

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Shufen Qin

    % Normalization
    PopObj = (PopObj - min(PopObj,[],1))./(max(PopObj,[],1)-min(PopObj,[],1));
    
    % Associate each solution to the closest reference vector
    Angle         = acos(1-pdist2(PopObj,V,'cosine'));
    [~,associate] = min(Angle,[],2);
    
    % Select the next population
    NV   = size(V,1);
    Next = zeros(1,NV);
    for i = unique(associate)'
        current  = find(associate==i);
        % Calculate the Euclidean distance
        de       = sqrt(sum(PopObj(current,:).^2,2));
        % Select the one with the minimum distance
        [~,best] = min(de);
        Next(i)  = current(best);
    end
    
    % Population for next generation
    index = Next(Next~=0);
end