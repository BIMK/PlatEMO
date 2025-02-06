function [Populations2,Masks2,Decs2,GV2,numK] = SubPopSimility(Populations,Masks,Decs,GV)
% Calculate the similarity between subpopulations

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Xiangyu Wang (email: xiangyu.wang@uni-bielefeld.de)

    K  = length(Populations);
    fs = zeros(1,K);
    for m  = 1 : K
        [~,fs(m)] = min(mean(Populations{m}.objs,2));
    end
    ss    = 0;
    flagmerge=zeros(K,1);
    for i = 1 : K-1
        if flagmerge(i,:)==0||flagmerge(i,:)==i
            for j = i+1 : K
                if flagmerge(j,:)==0
                        ss    = simility(Masks{i}(fs(i),:),Masks{j}(fs(j),:));
                        if ss>0.7
                            flagmerge(i,:)=i;
                            flagmerge(j,:)=i;
                        end
                else
                    continue;
                end
            end
        else
            continue;
        end
    end
    a    = find(flagmerge==0);
    b    = max(flagmerge);
    numK = size(a,1) + b;
    for k = unique(flagmerge)'
        temP2 = [];
        temD2 = [];
        temM2 = [];
        temG2 = [];
        current = find(flagmerge==k);
        if k ~= 0
            for h = 1 : size(current,1)
                temP  = Populations{current(h)};
                temP2 = [temP2,temP];
                temD  = Decs{current(h)};
                temD2 = [temD2;temD];
                temM  = Masks{current(h)};
                temM2 = [temM2;temM];
                temG  = GV{current(h)};
                temG2 = [temG2;temG];
            end
            Populations2{k+size(a,1)} = temP2;
            Decs2{k+size(a,1)}        = temD2;
            Masks2{k+size(a,1)}       = temM2;
            GV2{k+size(a,1)}          = temG2;
        else
            for h = 1 : size(current,1)
                temP = Populations{current(h)};
                temD = Decs{current(h)};
                temM = Masks{current(h)};
                temG = GV{current(h)};
                Populations2{h} = temP;
                Decs2{h} = temD;
                Masks2{h}= temM;
                GV2{h}   = temG;
            end
        end
    end
    Populations2(cellfun(@isempty,Populations2)) = [];
    Decs2(cellfun(@isempty,Decs2))   = [];
    Masks2(cellfun(@isempty,Masks2)) = [];
    GV2(cellfun(@isempty,GV2)) = [];
    numK = size(Populations2,2);
end

function s = simility(subPop1,subPop2)
	s = sum(subPop1&subPop2)/min(sum(subPop1),sum(subPop2));
end