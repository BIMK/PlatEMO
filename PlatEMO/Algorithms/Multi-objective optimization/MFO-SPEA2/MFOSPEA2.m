classdef MFOSPEA2 < ALGORITHM
% <2023> <multi> <real/integer/label/binary/permutation> <constrained>
% Multiform optimization framework based on SPEA2

%------------------------------- Reference --------------------------------
% R. Jiao, B. Xue, and M. Zhang. A multiform optimization framework for
% constrained multiobjective optimization. IEEE Transactions on
% Cybernetics, 2023, 53(8):5165-5177.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Ruwang Jiao

    methods
        function main(Algorithm, Problem)
            %% Generate random population for target task
            TargetPop  = Problem.Initialization();
            %% Generate random population for source task
            SourcePop  = TargetPop;   
            initialE   = max(max(0,TargetPop.cons), [], 1);
            initialE(initialE<1) = 1;
            [~, rankTP] = CalFitness(TargetPop.objs, TargetPop.cons);
            [~, rankSP] = CalFitness(SourcePop.objs, SourcePop.cons-initialE);

            %% Optimization
            while Algorithm.NotTerminated(TargetPop)
                epsn       = ReduceBoundary(initialE, ceil(Problem.FE/Problem.N), ceil(Problem.maxFE/Problem.N)-1);
                MatingPool = TournamentSelection(2, Problem.N, [rankTP', rankSP']);
                P          = [TargetPop, SourcePop];
                Offspring  = OperatorGA(Problem, P(MatingPool));
                %% Environmental selection for target task
                [TargetPop, rankTP] = EnvironmentalSelection([TargetPop, Offspring], Problem.N, [TargetPop.cons; Offspring.cons]);
                %% Environmental selection for source task
                [SourcePop, rankSP] = EnvironmentalSelection([SourcePop, Offspring], Problem.N, [SourcePop.cons-epsn; Offspring.cons-epsn]);
            end
        end
    end
end

function epsn = ReduceBoundary(eF, k, MaxK)
    %% Reduce the epsilon constraint boundary for source task
    z        = 1e-8;
    Nearzero = 1e-15;
    B        = MaxK./power(log((eF + z)./z), 1.0./10);
    B(B==0)  = B(B==0) + Nearzero;
    f        = eF.* exp( -(k./B).^10 );
    tmp      = find(abs(f-z) < Nearzero);
    f(tmp)   = f(tmp).*0 + z;
    epsn     = f - z;
    epsn(epsn<=0) = 0;
end 