function [ss,index] = SubPopSimility(Populations,Masks)
% Calculate the similarity between subpopulations

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    K  = length(Populations);
    fs = zeros(1,K);
    for m  = 1 : K
        [~,fs(m)] = min(mean(Populations{m}.objs,2));
    end
    ss    = 0;
    index = [];
    for i = 1 : K-1
        for j = i+1 : K
            if simility(Masks{i}(fs(i),:),Masks{j}(fs(j),:)) > ss
                ss    = simility(Masks{i}(fs(i),:),Masks{j}(fs(j),:));
                index = [i,j];
            end
        end
    end
end

function s = simility(subPop1,subPop2)
	s = sum(subPop1&subPop2)/min(sum(subPop1),sum(subPop2));
end