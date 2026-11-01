function [TrendOffspring,OrthogonalOffspring] = CollaborativeLearning(Problem,Swarm1,Swarm2,W,epsilon)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yubo Wang (email: 1202111230@cug.edu.cn)

    [TrendOffspring,Direction] = TrendLearning(Problem,Swarm1,W,epsilon);
    OrthogonalOffspring        = OrthogonalLearning(Problem,Swarm2,W,Direction,0);
end

function [SampleSolution,Direction] = TrendLearning(Problem,Population,RefV,epsilon)
    Lower = Problem.lower;
    Upper = Problem.upper;
    tau   = Problem.FE / Problem.maxFE;

    Nr              = size(RefV,1);
    LowerDirections = zeros(Nr, Problem.D);
    UpperDirections = zeros(Nr, Problem.D);
    CPopNum         = floor(Problem.N/Nr);
    Obj             = Population.objs;
    np              = size(Obj,1);
    Obj             = (Obj-repmat(min(Obj),np,1))./(repmat(max(Obj),np,1)-repmat(min(Obj),np,1));
    for i = 1 : Nr
        Cosine        = 1-pdist2(RefV(i,:),Obj,'cosine');
        [~,I]         = maxk(Cosine,CPopNum);
        Pop           = Population(I);
        PopCons       = Pop.cons;
        FeasibleIdx   = find(all(PopCons <= epsilon, 2));
        InfeasibleIdx = find(any(PopCons > epsilon, 2));
        if ~isempty(FeasibleIdx)
            FeasObjs      = Pop(FeasibleIdx).objs;
            [FrontNo, ~]  = NDSort(FeasObjs, 1);
            WinnerIndices = FeasibleIdx(FrontNo == 1);
            LoserIndices  = FeasibleIdx(FrontNo > 1);
            if isempty(LoserIndices)
                if ~isempty(InfeasibleIdx)
                    LoserIndices = InfeasibleIdx; 
                else
                    LoserIndices = WinnerIndices;
                end
            end
            if isempty(InfeasibleIdx)
                LoserIndices = WinnerIndices;
            end
        else
            [~, RankCons] = sort(sum(max(0, PopCons), 2));
            HalfPop       = floor(length(RankCons) / 2);
            WinnerIndices = RankCons(1 : HalfPop);        
            LoserIndices  = RankCons(HalfPop+1 : end); 
        end
        x_s = Pop(WinnerIndices(randi(length(WinnerIndices)))).decs;
        x_e = Pop(LoserIndices(randi(length(LoserIndices)))).decs;
        LowerGlobalDirection  = x_s - Lower;
        LocalDirection1       = x_s - x_e;
        UpperGlobalDirection  = x_s - Upper;
        LocalDirection2       = x_s - x_e;
        LowerDirections(i, :) = (1 - tau) * LowerGlobalDirection + tau * LocalDirection1;
        UpperDirections(i, :) = (1 - tau) *  UpperGlobalDirection + tau * LocalDirection2;
    end
    Directnorm = [sqrt(sum(LowerDirections.^2, 2));
                  sqrt(sum(UpperDirections.^2, 2))];
    Directnorm(Directnorm == 0) = 1e-8;
    Direction      = [LowerDirections; UpperDirections] ./ repmat(Directnorm, 1, Problem.D);
    Intervalmax    = sqrt(sum((Upper - Lower).^2, 2));
    Intervalmin    = 0;
    RandSample     = Intervalmin + rand(30, 2*Nr) * (Intervalmax - Intervalmin);
    SampleSolution = GenerateSampleSolution(Problem, RandSample, Direction, Nr);
end

function SampleSolution = OrthogonalLearning(Problem,Population,RefV,Direction,epsilon)
    Lower           = Problem.lower;
    Upper           = Problem.upper;
    Nr              = size(RefV,1);
    Intervalmax     = sqrt(sum((Upper - Lower).^2, 2));
    Intervalmin     = 0;
    RandSteps_Div_L = Intervalmin + rand(30, Nr) * (Intervalmax - Intervalmin);
    RandSteps_Div_U = Intervalmin + rand(30, Nr) * (Intervalmax - Intervalmin);
    [~, D]          = size(Population.decs);
    CPopNum         = floor(Problem.N/Nr);
    Obj             = Population.objs;
    np              = size(Obj,1);
    Obj             = (Obj-repmat(min(Obj),np,1))./(repmat(max(Obj),np,1)-repmat(min(Obj),np,1));
    SampleSolution  = [];
    for i = 1 : Nr
        Cosine      = 1-pdist2(RefV(i,:),Obj,'cosine');
        [~,I]       = maxk(Cosine,CPopNum);
        Pop         = Population(I);
        PopCons     = Pop.cons;
        FeasibleIdx = find(all(PopCons <= epsilon, 2));
        if ~isempty(FeasibleIdx)
            FeasObjs      = Pop(FeasibleIdx).objs;
            [FrontNo, ~]  = NDSort(FeasObjs, 1);
            WinnerIndices = FeasibleIdx(FrontNo == 1);
        else
            [~, RankCons] = sort(sum(max(0, PopCons), 2));
            WinnerIndices = RankCons(1); 
        end
        x                  = Pop(WinnerIndices(randi(length(WinnerIndices)))).decs;
        v_Conv_L_Norm      = Direction(i, :);
        v_Conv_U_Norm      = Direction(Nr + i, :);
        v_Ortho_L_Norm     = CreateOrthogonalVector(v_Conv_L_Norm, D);
        v_Ortho_U_Norm     = CreateOrthogonalVector(v_Conv_U_Norm, D);
        steps_L            = RandSteps_Div_L(:, i);
        steps_U            = RandSteps_Div_U(:, i);
        PopX_L             = repmat(x, 30, 1) + repmat(steps_L, 1, D) .* repmat(v_Ortho_L_Norm, 30, 1);
        PopX_U             = repmat(x, 30, 1) + repmat(steps_U, 1, D) .* repmat(v_Ortho_U_Norm, 30, 1);
        PopX               = [PopX_L; PopX_U];
        PopX               = max(min(PopX, Upper), Lower);
        SampleSolutiontemp = Problem.Evaluation(PopX);
        SampleSolution     = [SampleSolution, SampleSolutiontemp];
    end
end

function v_ortho_normalized = CreateOrthogonalVector(v_conv_normalized, D)
    z       = randn(1, D);
    v_proj  = (sum(z .* v_conv_normalized)) * v_conv_normalized;
    v_ortho = z - v_proj;

    norm_v_ortho = norm(v_ortho);
    if norm_v_ortho < 1e-8
        v_ortho      = randn(1, D); 
        v_proj       = (sum(v_ortho .* v_conv_normalized)) * v_conv_normalized;
        v_ortho      = v_ortho - v_proj;
        norm_v_ortho = norm(v_ortho) + 1e-8;
    end
    v_ortho_normalized = v_ortho / norm_v_ortho;
end

function SampleSolution = GenerateSampleSolution(Problem, RandSample, Direct, Nr)
    [NUM, ~]       = size(RandSample);
    SampleSolution = [];
    for i = 1 : NUM
        PopX_L = repmat(Problem.lower, Nr, 1) + repmat(RandSample(i, 1:Nr)', 1, Problem.D) .* Direct(1:Nr, :);
        PopX_U = repmat(Problem.upper, Nr, 1) + repmat(RandSample(i, Nr+1:end)', 1, Problem.D) .* Direct(Nr+1:end, :);
        PopX   = [PopX_L; PopX_U];
        PopX   = max(min(repmat(Problem.upper, size(PopX, 1), 1), PopX), repmat(Problem.lower, size(PopX, 1), 1));
        SampleSolutiontemp = Problem.Evaluation(PopX);
        SampleSolution     = [SampleSolution, SampleSolutiontemp];
    end
end