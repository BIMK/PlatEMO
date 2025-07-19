function [Population,FrontNo,CrowdDis] = AE_prediction(Problem,curr_NDS,his_NDS,curr_POS,NP)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    index1  = curr_NDS.decs;
    index2  = his_NDS.decs;
    [d, N1] = size(index1');
    [~, N2] = size(index2');
    if N1>N2
        [FrontNo,~] = NDSort(curr_NDS.objs,curr_NDS.cons,N1);
        CrowdDis    = CrowdingDistance(curr_NDS.objs,FrontNo);
        [~,index]   = sort(CrowdDis);
        curr_NDS    = curr_NDS(index(1:N2));
        NL = N2;
    else
        [FrontNo,~] = NDSort(his_NDS.objs,his_NDS.cons,N2);
        CrowdDis    = CrowdingDistance(his_NDS.objs,FrontNo);
        [~,index]   = sort(CrowdDis);
        b           = index(1:N1);
        his_NDS     = his_NDS(b);
        NL = N1;
    end
    index1 = curr_NDS.decs;
    index2 = his_NDS.decs;
    Q      = index2*index2';
    P      = index1*index2';
    lambda = 1e-5;
    reg    = lambda*eye(NL);
    reg(end,end) = 0;
    M    = P/(Q+reg);
    varM = M*index2;
    for i = 1 : NL
        for j = 1 : d
            var(i,j) = (index1(i,j)-varM(i,j)).^2;
        end
    end
    v = mean2(var);
    
    pre_solution = M*index1+v;
    POP          = Problem.Evaluation(pre_solution);
    curr_len     = length(POP);
    if curr_len > NP/2
        [FrontNo,~] = NDSort(POP.objs,POP.cons,N1);
        CrowdDis    = CrowdingDistance(POP.objs,FrontNo);
        [~,index]   = sort(CrowdDis);
        POP         = POP(index(1:NP/2));
    end
    Selected   = randperm(length(curr_POS),NP/2);
    POP2       = curr_POS(Selected);
    Population = [POP,POP2];
    if length(Population) < NP
        N    = NP-length(Population);
        POP3 = Problem.Initialization(N);
        Population = [Population,POP3];
    end
    [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Problem.N);
end