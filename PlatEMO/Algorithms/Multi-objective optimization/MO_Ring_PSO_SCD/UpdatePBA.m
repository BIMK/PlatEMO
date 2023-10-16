function PBA = UpdatePBA(Population,PBA,n_PBA)
% Update PBA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    for i = 1 : length(PBA)
        tempPBA = [PBA{i},Population(i)];
        [tempPBA,FrontNo,SpCrowdDis] = non_domination_scd_sort(tempPBA,n_PBA);
        [~,index] = sortrows([FrontNo;-SpCrowdDis]');
        PBA{i} = tempPBA(index);
    end
end