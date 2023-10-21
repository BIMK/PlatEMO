classdef L2SMEA < ALGORITHM
% <single> <real> <expensive>
% Linear subspace surrogate modeling assisted evolutionary algorithm
% NLinear --- 8 --- Number of linear subspaces

%------------------------------- Reference --------------------------------
% S. L, X. Zhang, Y. Tian, S. Yang, L. Zhang, and Y. Jin, Linear subspace
% surrogate modeling for large-scale expensive single/multi-objective
% optimization, IEEE Transactions on Evolutionary Computation, 2023.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm, Problem)
            %% Parameter setting
            NLinear = Algorithm.ParameterSet(8);

            %% Initialize an Archive
            BU       = Problem.upper;
            BD       = Problem.lower;
            D        = Problem.D;
            InitNum  = 2*D;
            Decs1    = lhsdesign(InitNum,D);
            Decs     = Decs1.*repmat(BU-BD,InitNum,1) + repmat(BD,InitNum,1);
            Archive  = Problem.Evaluation(Decs);
            TArchive = [Decs1,Archive.objs];

            %% Default Parameters of multi-objective optimizer
            PopSize = 40;
            MaxIter = 30;

            %% Default Parameters of cma-es optimizer
            XMean = rand(D,1);
            Sigma = 0.5;
            % Initialize dynamic (internal) strategy parameters and constants
            PathConv  = zeros(D,1);
            PathSig   = zeros(D,1);
            B         = eye(D,D);
            Diag      = ones(D,1);
            Conv      = B*diag(Diag.^2)*B';
            InvSqrtC  = B*diag(Diag.^-1)*B';
            EigenEval = 0;
            ChiD      = D^0.5*(1-1/(4*D)+1/21*D^2);
            % Selection
            Lambda  = 4*NLinear;
            Mu      = Lambda/2;
            Weights = log(Mu+1/2) - log(1:Mu)';
            Mu      = floor(Mu);
            Weights = Weights/sum(Weights);
            MuEff   = sum(Weights)^2/sum(Weights.^2);
            % Adaptation
            CumConv = (4+MuEff/D)/(D+4+2*MuEff/D);
            CumSig  = (MuEff+2)/(D+MuEff+5);
            CRate   = 2/((D+1.3)^2+MuEff);
            CMu     = min(1-CRate,2*(MuEff-2+1/MuEff)/((D+2)^2+MuEff));
            Damps   = 1 + 2*max(0,sqrt((MuEff-1)/(D+1))-1);

            %% Optimization
            while Algorithm.NotTerminated(Archive)
                %% Construct several linear space
                % Select reference points
                RefPoints = SelectRef(TArchive,NLinear,XMean,Sigma,B,Diag);
                % Construct process
                SubproblemList = ConstructLinearKRG(RefPoints,TArchive,NLinear);

                %% Optimize process
                % Optimize by multi-objective algorithm
                Candidates = MOEAOptimize(SubproblemList,PopSize,MaxIter);
                % Select candidate solutions
                Offspring = SelectCandidateKRG(SubproblemList,NLinear,D,Candidates);

                %% Update Archive
                [TArchive,Archive] = UpdateArchive(TArchive,Problem,Offspring, Archive);

                %% Update parameters in CMA-ES
                [~,Index] = sort(TArchive(:,D+1),'ascend');
                XOld      = XMean;
                XMean     = TArchive(Index(1:Mu),1:D)'*Weights;

                PathSig  = (1-CumSig)*PathSig + sqrt(CumSig*(2-CumSig))*MuEff*InvSqrtC*(XMean-XOld)/Sigma;
                HSig     = norm(PathSig)/sqrt(1-(1-CumSig)^(2*(Problem.FE-InitNum)/Lambda))/ChiD < 1.4+2/(D+1);
                PathConv = (1-CumConv)*PathConv + HSig*sqrt(CumConv*(2-CumConv)*MuEff)*(XMean-XOld)/Sigma;

                ArcTmp = (1/Sigma)*(TArchive(Index(1:Mu),1:D)'-repmat(XOld,1,Mu));
                Conv   = (1-CRate-CMu)*Conv + CRate*(PathConv*PathConv'+(1-HSig)*CumConv*(2-CumConv)*Conv) + CMu*ArcTmp*diag(Weights)*ArcTmp';

                Sigma = Sigma*exp((CumSig/Damps)*(norm(PathSig)/ChiD-1));
                if Problem.FE-InitNum-EigenEval > Lambda/(CRate+CMu)/D/10
                    EigenEval = Problem.FE-InitNum;
                    Conv      = triu(Conv) + triu(Conv,1)';
                    [B,Diag]  = eig(Conv);
                    Diag      = sqrt(diag(Diag));
                    InvSqrtC  = B*diag(Diag.^-1)*B';
                    if ~isreal(Diag)
                        Diag = abs(Diag);
                    end
                end
            end
        end
    end
end