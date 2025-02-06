function [SubPopulation,Dec,Mask,FrontNo,CrowdDis,NumTransUp,ALLTthetamid,MeanMaskSourceId]= ...
    Op_Transnodec(Problem,SubPopulation,Dec,Mask,FrontNo,CrowdDis,EachN,TF1,FitnessDec,pDecIndex,FitnessPop,SourceId,Buquan,NumTransUp,ALLTthetamid)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    TaskNum    = size(Problem.SubM,2);
    NSpop      = {};
    NSdec      = {};
    NSmask     = {};
    NSCrowdDis = {};
    MeanMask   = zeros(TaskNum,max(Problem.SubD));
    for i = 1 : TaskNum
        NSpop{i}      = SubPopulation{i}(FrontNo{i} == 1);
        NSdec{i}      = Dec{i}(FrontNo{i} == 1,:);
        NSdec{i}      = [NSdec{i},zeros(size(NSdec{i},1),Buquan(i))];
        NSmask{i}     = Mask{i}(FrontNo{i} == 1,:);
        NSmask{i}     = [NSmask{i},zeros(size(NSmask{i},1),Buquan(i))];
        NSCrowdDis{i} = CrowdDis{i}(FrontNo{i} == 1);
        MeanMask(i,:) = mean(NSdec{i}.*NSmask{i});
    end
    MeanMaskSourceId = [];
    
    OffDec  = {};
    OffMask = {};
    Pop     = {};
    sRatio  = zeros(1,TaskNum);
    CurrentTthetamid = zeros(1,TaskNum);
    for i = 1 : TaskNum
        [OffDec{i},OffMask{i},CurrentNumTransing(i)] = Operator(i,SourceId(i),Problem,NSpop{i},NSdec{i},NSmask{i},NSCrowdDis{i},TF1(i,:),FitnessDec{i},pDecIndex(i,:),FitnessPop(:,i),...
            NSpop{SourceId(i)},NSdec{SourceId(i)},NSmask{SourceId(i)},NSCrowdDis{SourceId(i)},TF1(SourceId(i),:),FitnessDec{SourceId(i)},pDecIndex(SourceId(i),:),FitnessPop(:,SourceId(i)),NumTransUp,ALLTthetamid);
        Skill    = i*ones(size(OffDec{i},1),1);
        Solution = [OffDec{i}.*OffMask{i},Skill];
        if size(Solution,1) > 0
            Pop{i} = Problem.Evaluation(Solution);
        else
            Pop{i} = [];
        end
        [NSpop{i},NSdec{i},NSmask{i},TFrontNo,TCrowdDis,sRatio(i)] = MOEAPSL_EnvironmentalSelection([NSpop{i},Pop{i}],[NSdec{i};OffDec{i}],[NSmask{i};OffMask{i}],length(NSpop{i}),NumTransUp(i));
        NumTransUp(i) = round(NumTransUp(i)*(1+sRatio(i))/2);
        Theta         = sum(NSmask{i}(find(TFrontNo ==1),:),2)';
        CurrentTthetamid(i) = mean(Theta);
    end
    ALLTthetamid = [ALLTthetamid;CurrentTthetamid];
    for i = 1 : TaskNum
        [SubPopulation{i},Dec{i},Mask{i},FrontNo{i},CrowdDis{i}] = SparseEA_EnvironmentalSelection([SubPopulation{i},Pop{i}],[Dec{i};OffDec{i}(:,1:Problem.SubD(i))],[Mask{i};OffMask{i}(:,1:Problem.SubD(i))],EachN);
    end
end

function [OffDec,OffMask,N] = Operator(TargetId,SourceId,Problem,Pop1,Dec1,Mask1,CrowdDis1,Fitness1,FitnessDec1,pDecIndex1,FitnessPop1,Pop2,Dec2,Mask2,CrowdDis2,Fitness2,FitnessDec2,pDecIndex2,FitnessPop2,NumTransUp,ALLTthetamid)
    [N,D] = deal(size(Dec1,1),Problem.SubD(TargetId));
    if N>NumTransUp(TargetId)
        N = NumTransUp(TargetId);
        TMatingPool = TournamentSelection(2,N,-CrowdDis1);
        TDec        = Dec1(TMatingPool,:);
        TMask       = Mask1(TMatingPool,:);
    else
        TDec  = Dec1;
        TMask = Mask1;
    end
    Tthetamid = ALLTthetamid(end,TargetId);

    SMatingPool = TournamentSelection(2,N,-CrowdDis2);
    SDec        = Dec2(SMatingPool,:);
    SMask       = Mask2(SMatingPool,:);
    
    S1index = sum(SMask)>0;

    OffMask = TMask;
    for i = 1 : N
        if sum(OffMask(i,:)) <= Tthetamid
            index0           = find(~TMask(i,:)&S1index);
            index0(index0>D) = [];
            minN             = floor(Tthetamid - sum(OffMask(i,:)));
            index            = index0(TStransMany(Fitness1(index0),Fitness2(index0),minN));
            OffMask(i,index) = 1;
        else
            index            = find(~TMask(i,:)&SMask(i,:));
            index(index>D)   = [];
            index            = index(TStrans(Fitness1(index),Fitness2(index)));
            OffMask(i,index) = SMask(i,index);
        end
    end

    for i = 1 : N
        index            = find(OffMask(i,:));
        index(index>D)   = [];
        index            = index(TS(-Fitness1(index)));
        OffMask(i,index) = 0;
    end

    if any(Problem.encoding~=4)
        OffDec = TDec;
        OffDec(:,Problem.encoding==4) = 1;
    else
        OffDec = ones(size(OffMask));
    end
end

function index = TStransMany(Fitness1,Fitness2,minN)
% Binary tournament selection

    if isempty(Fitness1)
        index = [];
    else
        index = TournamentSelection(2,minN,Fitness1,Fitness2);
    end
end

function index = TS(Fitness)
% Binary tournament selection

    if isempty(Fitness)
        index = [];
    else
        index = TournamentSelection(2,1,Fitness);
    end
end

function index = TStrans(Fitness1,Fitness2)
% Binary tournament selection

    if isempty(Fitness1)
        index = [];
    else
        index = TournamentSelection(2,1,Fitness1,Fitness2);
    end
end
