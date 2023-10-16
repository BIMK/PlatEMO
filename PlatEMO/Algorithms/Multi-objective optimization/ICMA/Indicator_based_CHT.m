function SelInd = Indicator_based_CHT(PopObj,IMatrix,W,N)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Jiawei Yuan

    % calculate the size of individuals in the promising areas
    IrFitness                         = min(IMatrix,[],2);
    Level1Index                       = find(IrFitness>=0);
    Len_Level1                        = length(Level1Index);

    if Len_Level1<=N
        [~,SortIndex]                 = sort(-IrFitness);
        SelInd                        = SortIndex(1:N);
    else
        % only focus on the solutions in the promising areas
        SelInd                        = Level1Index;
        PopObj                        = PopObj(Level1Index,:);
        IMatrix                       = IMatrix(Level1Index,Level1Index);

        [Num,M]                       = size(PopObj);
        NormW                         = W./repmat(sqrt(sum(W.^2,2)),1,M);
        NormPopObj                    = PopObj./repmat(sqrt(sum(PopObj.^2,2)),1,M);
        [~,ZoneIndex]                 = max(NormPopObj * NormW',[],2);
        Num_W                         = size(W,1);
        ZoneDensity                   = zeros(1,Num_W);
        zone.index                    = [];
        Zone                          = repmat(zone,1,Num_W);
        for j = 1:Num
            Zj                        = ZoneIndex(j);
            Zone(Zj).index            = [Zone(Zj).index,j];
            ZoneDensity(Zj)           = ZoneDensity(Zj) + 1;
        end

        [NDensity,SortIndex]          = sort(-ZoneDensity);
        Density                       = abs(NDensity);
        [Values,Neightboor]           = min(IMatrix,[],2);

        DelNum           = Num - N;
        Have_Delect      = zeros(1,DelNum);

        for i = 1:DelNum
            [MDen,MDInd] = max(Density);
            CandidateIndex           = Zone(SortIndex(MDInd)).index;

            [~,NowDel_Ind]           = min(Values(CandidateIndex));
            Del_Ind                  = CandidateIndex(NowDel_Ind);
            CandidateIndex(NowDel_Ind) = [];
            Have_Delect(i) = Del_Ind;
            IMatrix(Del_Ind,:) = Inf;
            IMatrix(:,Del_Ind) = Inf;
            Need_Updata = find(Neightboor==Del_Ind);
            L_Need=length(Need_Updata);
            if L_Need>0
                [Values(Need_Updata),Neightboor(Need_Updata)]=min(IMatrix(Need_Updata,:),[],2);
            end
            Values(Del_Ind) = Inf;

            Zone(SortIndex(MDInd)).index = CandidateIndex;
            Density(MDInd) = MDen - 1;

        end
        SelInd(Have_Delect) = [];
    end
end