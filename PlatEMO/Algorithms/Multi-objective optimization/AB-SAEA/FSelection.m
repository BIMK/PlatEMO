function index = FSelection(FunctionValue,V,theta0,Flag)
% The environmental selection of K-RVEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [N,M] = size(FunctionValue);
    NV    = size(V,1);
    
    %% Translate the population
    FunctionValue = FunctionValue - repmat(min(FunctionValue,[],1),N,1);
    
    %% Calculate the smallest angle value between each vector and others
    cosine = 1 - pdist2(V,V,'cosine');
    cosine(logical(eye(length(cosine)))) = 0;
    gamma = min(acos(cosine),[],2);

    %% Associate each solution to a reference vector
    Angle = acos(1-pdist2(FunctionValue,V,'cosine'));
    [~,associate] = min(Angle,[],2);
    
    %% Select one solution for each reference vector
    Next = zeros(1,NV);
    for i = unique(associate)'
        current = find(associate==i);
        if Flag==2
            % Calculate the APD value of each solution
             AD = (1+M*theta0*Angle(current,i)/gamma(i)).*sqrt(sum(FunctionValue(current,:).^2,2));
        else
             % Calculate the Angle value of each solution 
             AD = M*Angle(current,i)/gamma(i);   
        end
        % Select the one with the minimum AD value
        [~,best] = min(AD);
        Next(i)  = current(best);
    end
    % Population for next generation
    index = Next(Next~=0);
end