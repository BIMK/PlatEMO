function [SubPopulation,SubDec,SubMask,Rank] = InitRank(SubPopulation,SubDec,SubMask,Problem)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    Rank = {};
    SubPopulationCount = zeros(1,size(SubPopulation,2));
    for i = 1 : size(SubPopulation,2)
        SubPopulationCount(i) = size(SubPopulation{i},2);
    end
    for j = 1 : size(SubPopulation,2)
        [~,FrontNo{j},CrowdDis{j}] = EnSelection(SubPopulation{j},SubPopulationCount(j));
        [FrontNo{j},I1] = sort(FrontNo{j});
        CrowdDis{j} = CrowdDis{j}(I1);

        SubDec{j} = SubDec{j}(I1,:);

        SubMask{j} = SubMask{j}(I1,:);

        SubPopulation{j} = SubPopulation{j}(I1);
        count = zeros(1,max(FrontNo{j}));
        for i = 1 : max(FrontNo{j})
            for m = 1 : SubPopulationCount(j)
                if FrontNo{j}(m)==i
                    count(i) = count(i)+1;
                end
            end
        end
        NFindex = 1;
        for i = 1 : max(FrontNo{j})
            NFcount = count(i);
            [CrowdDis{j}(NFindex:NFindex+NFcount-1),I2] = sort(-CrowdDis{j}(NFindex:NFindex+NFcount-1));
            tempPop = SubPopulation{j}(NFindex:NFindex+NFcount-1);

            tempDec  = SubDec{j}(NFindex:NFindex+NFcount-1,:);
            tempMask = SubMask{j}(NFindex:NFindex+NFcount-1,:);

            tempPop = tempPop(I2);

            tempDec  = tempDec(I2,:);
            tempMask = tempMask(I2,:);

            SubPopulation{j}(NFindex:NFindex+NFcount-1) = tempPop;

            SubDec{j}(NFindex:NFindex+NFcount-1,:)  = tempDec;
            SubMask{j}(NFindex:NFindex+NFcount-1,:) = tempMask;

            NFindex = NFindex+NFcount;
        end
    end
    for i = 1 : size(SubPopulation,2)
        Rank{i} = 1:SubPopulationCount(i);
    end
end