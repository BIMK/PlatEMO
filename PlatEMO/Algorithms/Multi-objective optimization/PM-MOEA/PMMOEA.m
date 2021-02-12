classdef PMMOEA < ALGORITHM
% <multi> <real/binary> <large/none> <constrained/none> <sparse>
% Pattern mining based multi-objective evolutionary algorithm

%------------------------------- Reference --------------------------------
% Y. Tian, C. Lu, X. Zhang, F. Cheng, and Y. Jin, A pattern mining based
% evolutionary algorithm for large-scale sparse multi-objective
% optimization problems, IEEE Transactions on Cybernetics, 2020.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Population initialization
            if strcmp(Problem.encoding,'binary')
                Dec = ones(Problem.N,Problem.D);
            else
                Dec = unifrnd(repmat(Problem.lower,Problem.N,1),repmat(Problem.upper,Problem.N,1));
            end
            Mask = false(size(Dec));
            for i = 1 : Problem.N
                Mask(i,randperm(end,ceil(rand.^2*end))) = true;
            end
            Population = SOLUTION(Dec.*Mask);
            [Population,Dec,Mask,FrontNo] = EnvironmentalSelection(Population,Dec,Mask,Problem.N);

            %% Optimization
            MaxP = false(20,Problem.D);
            MinP = MaxP;
            while Algorithm.NotTerminated(Population)
                [MaxP,MinP,Nonzero] = POSMining(logical(Population(FrontNo==1).decs),MaxP,MinP,20);
                MatingPool = TournamentSelection(2,Problem.N,FrontNo);
                [OffDec,OffMask] = Operator(Dec(MatingPool,:),Mask(MatingPool,:),MaxP(:,Nonzero),MinP(:,Nonzero),Nonzero,Population.decs);
                if ~isempty(OffDec)
                    [Population,Dec,Mask,FrontNo] = EnvironmentalSelection([Population,SOLUTION(OffDec.*OffMask)],[Dec;OffDec],[Mask;OffMask],Problem.N);
                end
            end
        end
    end
end