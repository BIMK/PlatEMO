classdef MOL2SMEA < ALGORITHM
% <multi> <real> <expensive> <large/none>
% Multi-objective linear subspace surrogate modeling assisted evolutionary algorithm
% NLinear --- 8 --- Number of one-dimensional models

%------------------------------- Reference --------------------------------
% L. Si, X. Zhang, Y. Tian, S. Yang, L. Zhang, and Y. Jin, Linear subspace
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
        function main(Algorithm,Problem)
            %% Generate random population
            NLinear    = Algorithm.ParameterSet(8);
            NT         = 2*Problem.D;
            MaxFEs     = Problem.maxFE - NT;
            Dec        = lhsdesign(NT,Problem.D).*repmat(Problem.upper-Problem.lower,NT,1) + repmat(Problem.lower,NT,1);
            Population = Problem.Evaluation(Dec);
            Objs       = Population.objs;
            MaxObj     = max(Objs,[],1);
            MinObj     = min(Objs,[],1);
            TArchive   = [Dec,Population.objs];
            TArchive(:,1:Problem.D) = (TArchive(:,1:Problem.D)-repmat(Problem.lower,NT,1))./repmat(Problem.upper-Problem.lower,NT,1);
            
            %% Set parameters
            % Parameter of multi-objective algorithm
            PopSize = 40;
            MaxIter = 30;
            % Parameter of surrogate model
            CurrFEs = 0;
            
            %% Cluster number is 5
            CLTn     = 5;
            [W,CLTn] = UniformPoint(CLTn, Problem.M);
            % Parameter in CMA-ES
            XMeanCell    = cell(1, CLTn);
            SigmaCell    = cell(1, CLTn);
            PathConvCell = cell(1, CLTn);
            PathSigCell  = cell(1,CLTn);
            BCell        = cell(1,CLTn);
            DiagCell     = cell(1,CLTn);
            ConvCell     = cell(1,CLTn);
            InvSqrtCCell = cell(1,CLTn);
            NObjs        = (TArchive(:, end-Problem.M+1:end) - repmat(MinObj,size(TArchive,1),1))./repmat(MaxObj-MinObj,size(TArchive,1),1);
            [~,grp]      = min(pdist2(W, NObjs,'cosine'), [], 1);
            for k = 1 : CLTn
                XMeanCell{k}    = mean(TArchive(grp==k,1:Problem.D),1)';
                SigmaCell{k}    = 0.5;
                PathConvCell{k} = zeros(Problem.D,1);
                PathSigCell{k}  = zeros(Problem.D,1);
                BCell{k}        = eye(Problem.D,Problem.D);
                DiagCell{k}     = ones(Problem.D,1);
                ConvCell{k}     = BCell{k}*diag(DiagCell{k}.^2)*BCell{k}';
                InvSqrtCCell{k} = BCell{k}*diag(DiagCell{k}.^-1)*BCell{k}';
            end
            EigenEval = 0;
            ChiD      = Problem.D^0.5*(1-1/(4*Problem.D)+1/21*Problem.D^2);

            % Selection
            Lambda  = 4*NLinear;
            Mu      = Lambda/2;
            Weights = log(Mu+1/2) - log(1:Mu)';
            Mu      = floor(Mu);
            Weights = Weights/sum(Weights);
            MuEff   = sum(Weights)^2/sum(Weights.^2);

            % Adaptation
            CumConv = (4+MuEff/Problem.D)/(Problem.D+4+2*MuEff/Problem.D);
            CumSig  = (MuEff+2)/(Problem.D+MuEff+5);
            CRate   = 2/((Problem.D+1.3)^2+MuEff);
            CMu     = min(1-CRate,2*(MuEff-2+1/MuEff)/((Problem.D+2)^2+MuEff));
            Damps   = 1 + 2*max(0,sqrt((MuEff-1)/(Problem.D+1))-1);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                NObjs   = (TArchive(:, end-Problem.M+1:end) - repmat(MinObj,size(TArchive,1),1))./repmat(MaxObj-MinObj,size(TArchive,1),1);
                [~,grp] = min(pdist2(W, NObjs,'cosine'), [], 1);
                for k = 1 : CLTn
                    %% Construct several linear space
                    % Select reference points
                    RefPoints = SelectRef(NLinear,XMeanCell{k},SigmaCell{k},BCell{k},DiagCell{k});
                    tArc = TArchive(grp==k,:);
                    if sum(grp==k) < 30
                        disth = pdist2(W(k,:),NObjs,'cosine');
                        [~,indx] = sort(disth);
                        tArc = TArchive(indx(1:2*Problem.D),:);
                    end
                    % Construct process
                    SubproblemList = ConstructLinearKRGOS(RefPoints,tArc,NLinear, MaxObj,MinObj, W(k,:));
                    
                    %% Optimize process
                    % Optimize by multi-objective algorithm
                    Candidates = MOEAOptimizeOS(SubproblemList,PopSize,MaxIter);
                    % Select candidate solutions
                    Offspring = SelectCandidateKRGOS(SubproblemList,NLinear,Problem.D,Candidates);
                    Remain    = MaxFEs - CurrFEs;
                    
                    %% Update Archive
                    [TArchive,CurrFEs,tArc,Population, MaxObj,MinObj] = UpdateArchiveCLT(TArchive,Problem,Offspring,Remain,CurrFEs,tArc,Population,MaxObj,MinObj);
                    [FrontNo,~] = NDSort(tArc(:, Problem.D+1:end),size(tArc,1));
                    
                    if size(tArc,1) >= Mu
                        in = find(FrontNo==1);
                        id = find(FrontNo~=1);
                        if length(in) >= ceil(Mu/2)
                            if length(id) >= ceil(Mu/2)
                                sni     = randperm(length(in),ceil(Mu/2));
                                sdi     = randperm(length(id),ceil(Mu/2));
                                MeanDec = [tArc(in(sni),1:Problem.D);tArc(id(sdi),1:Problem.D)];
                            else
                                sni     = randperm(length(in),Mu-length(id));
                                MeanDec = [tArc(id,1:Problem.D);tArc(in(sni),1:Problem.D)];
                            end                    
                        else
                            sdi     = randperm(length(id),Mu-length(in));
                            MeanDec = [tArc(in,1:Problem.D);tArc(id(sdi),1:Problem.D)];
                        end
                    else
                        dist    = pdist2(W(k,:), NObjs,'cosine');
                        [~,cdt] = sort(dist,'ascend');
                        MeanDec = TArchive(cdt(1:Mu),1:Problem.D);
                    end

                    XOld         = XMeanCell{k};
                    XMeanCell{k} = MeanDec'*Weights;
                    
                    PathSigCell{k}  = (1-CumSig)*PathSigCell{k} + sqrt(CumSig*(2-CumSig))*MuEff*InvSqrtCCell{k}*(XMeanCell{k}-XOld)/SigmaCell{k};
                    HSig            = norm(PathSigCell{k})/sqrt(1-(1-CumSig)^(2*CurrFEs/Lambda))/ChiD < 1.4+2/(Problem.D+1);
                    PathConvCell{k} = (1-CumConv)*PathConvCell{k} + HSig*sqrt(CumConv*(2-CumConv)*MuEff)*(XMeanCell{k}-XOld)/SigmaCell{k};
                    ArcTmp          = (1/SigmaCell{k})*(MeanDec'-repmat(XOld,1,Mu));
                    
                    ConvCell{k} = (1-CRate-CMu)*ConvCell{k} + CRate*(PathConvCell{k}*PathConvCell{k}'+(1-HSig)*CumConv*(2-CumConv)*ConvCell{k}) + CMu*ArcTmp*diag(Weights)*ArcTmp';
                    
                    SigmaCell{k} = SigmaCell{k}*exp((CumSig/Damps)*(norm(PathSigCell{k})/ChiD-1));
                    if CurrFEs - EigenEval > Lambda/(CRate+CMu)/Problem.D/10
                        EigenEval       = CurrFEs;
                        ConvCell{k}     = triu(ConvCell{k}) + triu(ConvCell{k},1)';
                        [BCell{k},DiagCell{k}] = eig(ConvCell{k});
                        DiagCell{k}     = sqrt(diag(DiagCell{k}));
                        InvSqrtCCell{k} = BCell{k}*diag(DiagCell{k}.^-1)*BCell{k}';
                        if ~isreal(DiagCell{k})
                            DiagCell{k} = real(DiagCell{k});
                        end
                    end
                end
            end
        end
    end
end