function WHVLoss = CalWHVLoss(PopObj,FrontNo,wz,AA,RA)
% Calculate the weighted hypervolume (WHV) loss of each solution

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Calculate the weight of each solution for WHV calculation
    Weight = ones(1,size(PopObj,1));
    if nargin > 2 && ~isempty(wz)
        for i = 1 : size(PopObj,1)
            if any(all(repmat(PopObj(i,:),size(AA,1),1)<=AA,2))
                % The solution belongs to the preferred region (Pr)
                Weight(i) = wz(3);
            elseif any(all(repmat(PopObj(i,:),size(RA,1),1)>RA,2))
                % The solution belongs to the dominated region (Do)
                Weight(i) = wz(1);
            else
                % The solution belongs to the no preference information
                % region (In)
                Weight(i) = wz(2);
            end
        end
    end
    
    %% Calculate the WHV loss of each solution front by front
    WHVLoss  = zeros(1,size(PopObj,1));
    RefPoint = max(PopObj,[],1) + 0.1;
    for f = setdiff(unique(FrontNo),inf)
        current  = find(FrontNo==f);
        totalWHV = CalWHV(PopObj(current,:),RefPoint,Weight(current));
        for i = 1 : length(current)
            drawnow();
            currenti           = current([1:i-1,i+1:end]);
            WHVLoss(current(i))= totalWHV - CalWHV(PopObj(currenti,:),RefPoint,Weight(currenti));
        end
    end
end