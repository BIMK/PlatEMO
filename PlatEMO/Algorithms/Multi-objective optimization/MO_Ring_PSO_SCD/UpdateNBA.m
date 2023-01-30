function NBA = UpdateNBA(NBA,n_NBA,PBA)
% Update NBA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    for i = 1:length(NBA)
        if i == 1
            tempNBA = [PBA{end},PBA{1},PBA{2}];
        elseif i == length(NBA)
            tempNBA = [PBA{end-1},PBA{end},PBA{1}];
        else
            tempNBA = [PBA{i-1},PBA{i},PBA{i+1}];
        end
        tempNBA = [tempNBA,NBA{i}];
        [tempNBA,FrontNo,SpCrowdDis] = non_domination_scd_sort(tempNBA,n_NBA);
        [~,index] = sortrows([FrontNo;-SpCrowdDis]');
        NBA{i} = tempNBA(index);
    end
end