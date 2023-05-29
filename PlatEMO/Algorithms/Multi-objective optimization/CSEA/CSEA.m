classdef CSEA < ALGORITHM
% <multi/many> <real/integer> <expensive>
% Classification based surrogate-assisted evolutionary algorithm
% k    ---    6 --- Number of reference solutions
% gmax --- 3000 --- Number of solutions evaluated by surrogate model

%------------------------------- Reference --------------------------------
% L. Pan, C. He, Y. Tian, H. Wang, X. Zhang, and Y. Jin, A classification
% based surrogate-assisted evolutionary algorithm for expensive
% many-objective optimization, IEEE Transactions on Evolutionary
% Computation, 2019, 23(1): 74-88.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [k,gmax] = Algorithm.ParameterSet(6,3000);

            %% Initalize the population by Latin hypercube sampling
            N          = min(11*Problem.D-1,109);
            PopDec     = UniformPoint(N,Problem.D,'Latin');
            Population = Problem.Evaluation(repmat(Problem.upper-Problem.lower,N,1).*PopDec+repmat(Problem.lower,N,1));
            Arc        = Population;
            
            %% Initialize the network
            hiddenLayerSize = ceil(Problem.D*2);
            layers = [featureInputLayer(Problem.D,'Normalization', 'zscore')
                    fullyConnectedLayer(hiddenLayerSize)
                    batchNormalizationLayer
                    reluLayer
                    fullyConnectedLayer(1)
                    sigmoidLayer
                    regressionLayer];

            maxEpochs = 400;
            miniBatchSize = 32;
            options = trainingOptions('adam', ...
                        'ExecutionEnvironment','auto', ...
                        'MaxEpochs',maxEpochs, ...
                        'MiniBatchSize',miniBatchSize, ...
                        'Shuffle','every-epoch', ...
                        'Plots','none', ...
                        'Verbose',false);%, ... Plots','none'

            %% Optimization
            while Algorithm.NotTerminated(Arc)
                % Select reference solutions and preprocess the data
                Ref    = RefSelect(Population,k);
                Input  = Population.decs;  
                Output = GetOutput(Population.objs,Ref.objs); 
                rr     = sum(Output)/length(Output);
                tr     = min(rr,1-rr)*0.5;
                [TrainIn,TrainOut,TestIn,TestOut] = DataProcess(Input,Output);
                net = trainNetwork(TrainIn,TrainOut-0,layers,options);

                % Error rates calculation
                TestPre = predict(net,TestIn);
                IndexGood = TestOut==1;
                p0 = sum(abs((TestOut(IndexGood)-TestPre(IndexGood))))/sum(IndexGood);
                p1 = sum(abs((TestOut(~IndexGood)-TestPre(~IndexGood))))/sum(~IndexGood);

                % Surrogate-assisted selection and update the population
                Next = SurrogateAssistedSelection(Problem,net,p0,p1,Ref,Population.decs,gmax,tr);
                if ~isempty(Next)
                    Arc = [Arc,Problem.Evaluation(Next)];
                end
                Population = RefSelect(Arc,Problem.N);
            end
        end
    end
end