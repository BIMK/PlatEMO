function MPAES(Global)
% <algorithm> <H-N>
% M-PAES: A Memetic Algorithm for Multiobjective Optimization
% l_fails   ---  5 --- Maximum number of consecutive failing local moves
% l_opt     --- 10 --- Maximum number of local moves
% cr_trials --- 20 --- Number of crossover trials
% div       --- 10 --- The number of divisions in each objective
% operator         --- FEP

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [l_fails,l_opt,cr_trials,div] = Global.ParameterSet(5,10,20,10);
    
    %% Generate random population
    P = Global.Initialization();
    G = P(NDSort(P.objs,1)==1);
    
    %% Optimization
    while Global.NotTermination(G)
        for i = 1 : Global.N
            H = G(~all(G.objs<=repmat(P(i).obj,length(G),1),2));
            [P(i),G] = PAES(Global,P(i),G,[H,P(i)],l_fails,l_opt,div);
        end
        P1(1:Global.N) = INDIVIDUAL();
        for i = 1 : Global.N
            for r = 1 : cr_trials
                Combine = [P,G];
                parents = Combine(randperm(length(Combine),2));
                c       = Global.Variation(parents,1,@EAreal,{1,20,0,0});
                [G,dominated,GCrowd,cCrowd,pCrowd] = UpdateArchive(G,c,parents,Global.N,div);
                if ~dominated && any(cCrowd<=pCrowd)
                    break;
                end
            end
            if dominated
                c = G(TournamentSelection(2,1,GCrowd));
            end
            P1(i) = c;
        end
        P = P1;
    end
end